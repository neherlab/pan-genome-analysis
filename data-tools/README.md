# panx-data: rebuild the published panX data of 2022

The viewer at `pangenome.org` loads its data from the S3 bucket `pangenome.org-data-master`: one folder per dataset, every file stored plain and as a `.gz` copy, 23.6 million objects for the 264 datasets of the 2022 pipeline run (PanXX, a Nextflow rewrite of this pipeline). The pipeline output of that run is archived on sciCORE, the cluster of the University of Basel, as one tarball per dataset. `panx-data` rebuilds the bucket's 2022 datasets from these archives, byte for byte, and checks the result against listings of the bucket. A download of the whole bucket is not needed (it would cost about $26: 40.8 million requests and 101.5 GB of transfer). Only 432 objects that S3 changed after the upload come from S3, and with the listings the steps below cost about $0.01 in requests plus at most $0.04 in transfer.

The tool is separate from the Python 2.7 pipeline in the parent folder and needs only the Python standard library.

## Requirements

- Python 3.11 or later, `jq`, `rsync`, and the AWS CLI
- `pigz` 2.8, linked with zlib 1.3: the `.gz` copies in the bucket are `pigz` output, and only these versions produce the same bytes. `rebuild` checks the versions before it starts
- A sciCORE account with read access to `/scicore/home/neher/GROUP/data/panXX/`, from the University network or through its VPN
- AWS credentials of the lab account with read access to `pangenome.org-data-master`. The steps below cost about $0.01 in requests and at most $0.04 in transfer ($0.09 per GB once the free 100 GB of the month are used up)
- About 210 GB of free disk space and 24 million free inodes

## Snapshot folder

`rebuild` and `verify` work on a snapshot folder. `src/` holds the inputs, each at the path that it has at its source, and stays unchanged. `rebuild` writes `out/`.

```
<snapshot>/
|-- src/
|   |-- MD5SUMS                                      checksum of every file in src/
|   |-- scicore/scicore/home/neher/GROUP/data/panXX/ copy of /scicore/home/neher/GROUP/data/panXX/ on sciCORE (no AWS cost):
|   |   |-- 13-07-2022/<dataset>.tar.gz              July 2022 run: 273 tarballs, 7.8 GB
|   |   |-- 04-11-2022/<dataset>.tar.gz              November 2022 rerun of 52 datasets, 2.3 GB
|   |   |-- 04-11-2022/zip_archives_update.zip       web/ files and download zips of the rerun
|   |   `-- zip_13-07-2022/<dataset>/*.zip           download zips of the July run, 4.6 GB
|   |-- s3/pangenome.org-data-master/<key>           objects copied from the bucket, at their keys (432 objects, 1 MB, about $0.0002)
|   `-- s3-listings/<date>/pangenome.org-data-master/  listings of the bucket (about $0.008 in requests, about 0.4 GB, $0 to $0.04 in transfer); panx_data/listings.py describes the files
`-- out/
    |-- pangenome.org-data-master/                   the rebuilt bucket
    |-- manifest.tsv.gz                              every file: path, size, MD5, source (scicore or s3)
    |-- deviations.tsv                               every difference from the listings, with its class
    |-- build-summary.tsv                            files, bytes, and time per dataset
    `-- BUILD_INFO.json                              tool commit, tool versions, counts, time
```

## Manual steps: from the sources to the checked rebuild

Run the steps in `bash` from this folder (`data-tools/`). Set the variables first:

```bash
snapshot=<snapshot folder>
scicore_user=<sciCORE user name>
bucket=pangenome.org-data-master
archives="${snapshot}/src/scicore/scicore/home/neher/GROUP/data/panXX"
objects="${snapshot}/src/s3/${bucket}"
listings="${snapshot}/src/s3-listings/$(date +%F)/${bucket}"
```

### 1. Copy the archives from sciCORE

Copy the three run folders through the transfer node, then compare names and sizes with sciCORE. Expect 874 files and 14.7 GB. This step makes no AWS requests and has no AWS cost; the 14.7 GB come over the University network.

```bash
mkdir -p "${archives}"
for run in 13-07-2022 04-11-2022 zip_13-07-2022; do
  rsync -a "${scicore_user}@transfer12.scicore.unibas.ch:/scicore/home/neher/GROUP/data/panXX/${run}" "${archives}/"
done
ssh "${scicore_user}@login12.scicore.unibas.ch" \
  'cd /scicore/home/neher/GROUP/data/panXX && find 13-07-2022 04-11-2022 zip_13-07-2022 -type f -printf "%p\t%s\n" | sort' > scicore-files.tsv
(cd "${archives}" && find 13-07-2022 04-11-2022 zip_13-07-2022 -type f -printf '%p\t%s\n' | sort) | diff scicore-files.tsv - && echo "archives equal sciCORE"
```

### 2. List the bucket

The rebuild takes the list of datasets from these listings, and the check compares against them. 273 LIST requests, about $0.0015, and about 2 MB of responses, below $0.001 in transfer:

```bash
mkdir -p "${listings}/top-level" "${listings}/full"
list() { aws s3api list-objects-v2 --bucket "${bucket}" --page-size 1000 "$@"; }
list --delimiter / --output json > "${listings}/root.json"
list --prefix web/ --output json > "${listings}/web.json"
list --prefix dataset/ --delimiter / --output json > "${listings}/dataset.json"
for dataset in $(jq -r '.CommonPrefixes[].Prefix | split("/")[1]' "${listings}/dataset.json"); do
  list --prefix "dataset/${dataset}/" --delimiter / --output json > "${listings}/top-level/${dataset}.json"
done
```

The datasets that a second August pass completed on 2022-08-26 need a full listing, because only the ETags tell which `.gz` bytes the bucket holds there (17 datasets, 1,191 LIST requests, about $0.0064, plus about 0.4 GB of responses, which cost up to $0.04 in transfer when the free 100 GB of the month are used up):

```bash
for dataset in $(jq -r '.Contents[]? | select(.LastModified | startswith("2022-08-26")) | .Key | split("/")[1]' "${listings}"/top-level/*.json | sort -u); do
  [[ -f "${archives}/13-07-2022/${dataset}.tar.gz" ]] || continue  # the older datasets are not rebuilt
  list --prefix "dataset/${dataset}/" --query 'Contents[].[Key,Size,ETag,LastModified]' --output text \
    > "${listings}/full/${dataset}.key-size-etag-modified.tsv"
done
```

### 3. Copy the objects that the archives cannot produce

432 GET requests, about $0.0002, and 1 MB of transfer, below $0.001:

- **The bucket root**: `index.json`, `not_found.json`, and their `.gz` copies
- **Later changes**: every top-level object of 2022-09-09, which are the `metainfo.tsv` files and their `.gz` copies of 212 datasets. S3 replaced them with a new column layout after the upload, and the tarballs hold the older version
- **4 `.gz` copies with another gzip header time**: the same content, but other bytes than the rebuild produces

```bash
{
  printf '%s\n' index.json index.json.gz not_found.json not_found.json.gz
  jq -r '.Contents[]? | select(.LastModified | startswith("2022-09-09")) | .Key' "${listings}"/top-level/*.json
  printf '%s\n' \
    dataset/Bacillus_velezensis/geneCluster/d66a12ef3a081038_refined_na_aln.fa.gz \
    dataset/Citrobacter_portucalensis/geneCluster/9671524fd5b514c4_refined_aa_aln.fa.gz \
    dataset/Citrobacter_portucalensis/geneCluster/fba6d8c25065b05f_refined_tree.json.gz \
    dataset/Phaeobacter_inhibens/geneCluster/be2a1c40c0196838_refined_aa_aln_reduced.fa.gz
} > s3-keys.txt
python3 -m panx_data fetch-s3 --keys s3-keys.txt --dest "${objects}"
```

### 4. Write the checksums and protect `src/`

```bash
(cd "${snapshot}/src" && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum) > "${snapshot}/MD5SUMS.tmp"
mv "${snapshot}/MD5SUMS.tmp" "${snapshot}/src/MD5SUMS"
chmod -R a-w "${snapshot}/src"
```

### 5. Rebuild and check

The rebuild and the check make no AWS requests. The rebuild takes about 8 minutes on 20 cores and writes 23,568,036 files (120.6 GB, about 190 GB on disk):

```bash
python3 -m panx_data rebuild "${snapshot}" --jobs 16
python3 -m panx_data verify "${snapshot}"
```

`verify` must report only known classes. For the bucket as listed on 2026-11-11, expect these counts (from the check of that day, limited to the listings of step 2):

| Class                            |  Count | Meaning                                                            |
| -------------------------------- | -----: | ------------------------------------------------------------------ |
| `missing-in-s3`                  | 16,053 | the upload left the object out; `out/` has the correct file        |
| `s3-empty`, `s3-shorter`         | 5 + 10 | the upload broke the object; `out/` has the correct file           |
| `only-in-s3:stale-archive`       |    434 | an old `*_gene_alignments.tar.gz` that nothing links               |
| `only-in-s3:dataset-not-rebuilt` |     74 | top-level objects of the 6 older datasets, which the archives lack |

A `differs` row means that S3 holds other bytes for a key. Copy that key into `src/s3/` as in step 3 (one GET request per key, $0.43 per 1M, plus its size in transfer), write `MD5SUMS` again, delete `out/`, and rebuild.

## Commands

```bash
python3 -m panx_data rebuild <snapshot>           # build <snapshot>/out/; --jobs N, --dataset NAME, --out DIR
python3 -m panx_data verify <snapshot>            # write out/deviations.tsv; exit 1 on an unknown difference
python3 -m panx_data fetch-s3 --keys FILE --dest <snapshot>/src/s3/pangenome.org-data-master
```

`rebuild` resumes an interrupted run and refuses to overwrite a complete one. `fetch-s3` uses the AWS CLI with the caller's credentials. It costs one GET request per object ($0.43 per 1M in `eu-central-1`) plus the transfer of the objects ($0.09 per GB above the free 100 GB per month).

## How the bucket was filled

`panx_data/rebuild.py` describes the upload of 2022 that the rebuild repeats: the decompression and `.gz` copy of every file, the alignment zips, the August rule that kept the pipeline's original `.gz` file when the new copy had the same size, and the objects that S3 changed later.

`verify` sorts every difference into a class. The known classes describe the bucket, not the rebuild:

- `missing-in-s3`: the upload left the object out
- `s3-empty`, `s3-shorter`: the upload broke the object
- `only-in-s3:stale-archive`: an old `*_gene_alignments.tar.gz` that nothing links
- `only-in-s3:dataset-not-rebuilt`: a dataset that the archives lack

Any other class means that the rebuild does not reproduce the bucket.

## Tests

```bash
uv run --group dev pytest
```
