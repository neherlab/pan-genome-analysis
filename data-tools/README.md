# panx-data: from the panX pipeline archives to the data of pangenome.org

The viewer at `pangenome.org` reads its species data from a data host. Two layouts exist:

- **The old layout**: the S3 bucket `pangenome.org-data-master`, one folder per dataset, every file stored plain and as a `.gz` copy, 23.6 million objects for the 264 datasets of the 2022 pipeline run (PanXX, a Nextflow rewrite of this pipeline)
- **Data format**: one immutable release folder per data version, 8 objects per dataset and about 2,200 objects per release, served at `data.pangenome.org`. `docs/data-format.md` of the app repository `neherlab/pan-genome-visualization` defines the format

`panx-data` builds both from the same inputs: the pipeline output of the 2022 run and the 2018 archive of six older datasets, both archived on sciCORE, the cluster of the University of Basel, plus a few objects of the old bucket. A download of the whole bucket is not needed (it would cost about $26: 40.8 million requests and 101.5 GB of transfer).

The tool is separate from the Python 2.7 pipeline in the parent folder and needs only the Python standard library.

## The path from the sources to the running app

1. Copy the archives from sciCORE and list the old bucket (steps 1 to 3 below)
2. Write the checksums of the inputs (step 4)
3. Check the inputs against the old bucket with `rebuild` and `verify` (step 5). A release does not need this check again
4. Build and check a release (step 6)
5. Check the release with the format checker of the app, and run the app on it (step 7)
6. Publish the release to the data host (step 8)

## Requirements

- [mise](https://mise.jdx.dev) and Docker. `.config/mise.toml` pins Python (a python-build-standalone build with its own zlib), `just`, `uv`, `ruff`, the AWS CLI, and the lints. `dev/docker/run` runs a command in the development container: Debian with these tools and `pigz` 2.8 with zlib 1.3, which the `.gz` copies of the old bucket need. The container has the same Python and zlib as `mise install` on the host, so a release has the same bytes on both
- A sciCORE account with read access to `/scicore/home/neher/GROUP/data/`, from the University network or through its VPN
- AWS credentials: read access to `pangenome.org-data-master` for steps 2 and 3, and the deploy role of the data host for step 8. The AWS prices below are those of `eu-central-1`: LIST and PUT $5.40, GET and HEAD $0.43 per 1M requests, transfer to the internet free for the first 100 GB of the month and $0.09 per GB above, storage in S3 Standard $0.0245 per GB-month
- About 230 GB of free disk space and 24 million free inodes for the old layout; a release needs about 14 GB more

Run every command from this folder (`data-tools/`), after `mise trust .config/mise.toml && mise install`. `dev/docker/run` mounts the snapshot folder that `PANX_SNAPSHOT` names, with `src/` read-only and `out/` writable, and a folder that `PANX_OUT` names, writable. It forwards the `AWS_*` variables.

## Snapshot folder

All commands work on a snapshot folder. `src/` holds the inputs, each at the path that it has at its source, and stays unchanged. The commands write `out/`.

```
<snapshot>/
|-- src/
|   |-- MD5SUMS                                        checksum of every file in src/
|   |-- scicore/scicore/home/neher/GROUP/data/
|   |   |-- panXX/                                     the 2022 run (no AWS cost):
|   |   |   |-- 13-07-2022/<dataset>.tar.gz            July 2022 run: 273 tarballs, 7.8 GB
|   |   |   |-- 04-11-2022/<dataset>.tar.gz            November 2022 rerun of 52 datasets, 2.3 GB
|   |   |   |-- 04-11-2022/zip_archives_update.zip     web/ files, the dataset list, and download zips of the rerun
|   |   |   `-- zip_13-07-2022/<dataset>/*.zip         download zips of the July run, 4.6 GB
|   |   `-- panX_archive/leftover_data/<name>.tar.gz   2018 archive of five older datasets, 3.5 GB (no AWS cost)
|   |-- s3/pangenome.org-data-master/<key>             objects copied from the bucket, at their keys (about 2 GB)
|   `-- s3-listings/<date>/pangenome.org-data-master/  listings of the bucket; panx_data/listings.py describes the files:
|       |-- root.json, web.json, dataset.json          the bucket root, web/, and the dataset folders
|       |-- top-level/<dataset>.json                   the top-level objects of each dataset
|       |-- full/<dataset>.*.tsv[.gz]                  every object of the datasets that need it (step 2)
|       `-- gene-cluster/<dataset>.*.tsv.gz            the geneCluster/ objects of the five older datasets of the 2018 archive
`-- out/
    |-- pangenome.org-data-master/                     the rebuilt old bucket, written by `rebuild`
    |-- manifest.tsv.gz                                every file of the old layout: path, size, MD5, source (scicore or s3)
    |-- deviations.tsv                                 every difference from the listings, with its class
    |-- build-summary.tsv                              files, bytes, and time per dataset
    |-- BUILD_INFO.json                                tool commit, tool versions, counts, time
    |-- data.pangenome.org/                            data root, written by `release-build`:
    |   |-- releases.json                              the releases of the data root and the default release
    |   `-- <release id>/                              index.json, manifest.json, and one folder per dataset
    `-- release-check.txt                              the findings of `release-check`
```

## Steps

Set the variables first:

```bash
snapshot=<snapshot folder>
scicore_user=<sciCORE user name>
scicore_transfer=<sciCORE transfer node>
scicore_login=<sciCORE login node>
bucket=pangenome.org-data-master
data=/scicore/home/neher/GROUP/data
archives="${snapshot}/src/scicore/scicore/home/neher/GROUP/data"
objects="${snapshot}/src/s3/${bucket}"
listings="${snapshot}/src/s3-listings/$(date +%F)/${bucket}"
export PANX_SNAPSHOT="${snapshot}"
```

### 1. Copy the archives from sciCORE

Copy the three folders of the 2022 run and the five tarballs of the 2018 archive through the transfer node, then compare names and sizes with sciCORE. Expect 879 files and 18.2 GB. This step makes no AWS requests and has no AWS cost.

```bash
mkdir -p "${archives}/panXX" "${archives}/panX_archive/leftover_data"
for run in 13-07-2022 04-11-2022 zip_13-07-2022; do
  rsync -a "${scicore_user}@${scicore_transfer}:${data}/panXX/${run}" "${archives}/panXX/"
done
for name in Enterobacteriales Pmnawg Pseudomonadales S_pneumoniae616 Vibrionales; do
  rsync -a "${scicore_user}@${scicore_transfer}:${data}/panX_archive/leftover_data/${name}.tar.gz" "${archives}/panX_archive/leftover_data/"
done
files() {
  { find panXX/13-07-2022 panXX/04-11-2022 panXX/zip_13-07-2022 -type f -printf '%p\t%s\n'
    find panX_archive/leftover_data -maxdepth 1 -type f -name '*.tar.gz' -printf '%p\t%s\n'; } | sort
}
ssh "${scicore_user}@${scicore_login}" "cd ${data} && $(declare -f files) && files" > scicore-files.tsv
(cd "${archives}" && files) | diff scicore-files.tsv - && echo "archives equal sciCORE"
```

`Pmnawg.tar.gz` is the dataset `P_marinus_meta`.

### 2. List the bucket

The rebuild takes the list of datasets from these listings, and the check compares against them. 273 LIST requests, about $0.0015, and about 2 MB of responses, below $0.001 in transfer:

```bash
mkdir -p "${listings}/top-level" "${listings}/full" "${listings}/gene-cluster"
list() { aws s3api list-objects-v2 --bucket "${bucket}" --page-size 1000 "$@"; }
list --delimiter / --output json > "${listings}/root.json"
list --prefix web/ --output json > "${listings}/web.json"
list --prefix dataset/ --delimiter / --output json > "${listings}/dataset.json"
for dataset in $(jq -r '.CommonPrefixes[].Prefix | split("/")[1]' "${listings}/dataset.json"); do
  list --prefix "dataset/${dataset}/" --delimiter / --output json > "${listings}/top-level/${dataset}.json"
done
```

The datasets that a second August pass completed on 2022-08-26 need a full listing, because only the ETags tell which `.gz` bytes the bucket holds there (17 datasets, 1,191 LIST requests, about $0.0064, plus about 0.4 GB of responses, which cost up to $0.04 in transfer):

```bash
for dataset in $(jq -r '.Contents[]? | select(.LastModified | startswith("2022-08-26")) | .Key | split("/")[1]' "${listings}"/top-level/*.json | sort -u); do
  [[ -f "${archives}/panXX/13-07-2022/${dataset}.tar.gz" ]] || continue  # the older datasets come from step 3
  list --prefix "dataset/${dataset}/" --query 'Contents[].[Key,Size,ETag,LastModified]' --output text \
    > "${listings}/full/${dataset}.key-size-etag-modified.tsv"
done
```

The five older datasets of the 2018 archive need a listing of their `geneCluster/` folders, which compares each file of the archive with the served object (1.2 million keys, 1,222 LIST requests, about $0.0066, plus about 0.3 GB of responses, up to $0.03 in transfer):

```bash
for dataset in S_pneumoniae616 P_marinus_meta Enterobacteriales Pseudomonadales Vibrionales; do
  list --prefix "dataset/${dataset}/geneCluster/" --max-items 2000000 \
    --query 'Contents[].[Key,Size,ETag,LastModified]' --output text |
    gzip -n > "${listings}/gene-cluster/${dataset}.key-size-etag-modified.tsv.gz"
done
```

### 3. Copy the objects that the archives cannot produce

- **The bucket root**: `index.json`, `not_found.json`, and their `.gz` copies
- **Later changes of the 2022 datasets**: every top-level object of 2022-09-09, which are the `metainfo.tsv` files and their `.gz` copies of 212 datasets. S3 replaced them with a new column layout after the upload, and the tarballs hold the older version
- **4 `.gz` copies with another gzip header time**: the same content, but other bytes than the rebuild produces
- **`Neisseria_gonorrhoeae_ECDC`**: an older dataset that no archive holds, copied whole (32,020 objects)
- **Newer files of the older datasets of the 2018 archive**: the files whose S3 content differs from the archive and is newer, by the listings of step 2 (MD5 against ETag): the 10 alignment archives `all_gene_alignments.tar.gz` and `core_gene_alignments.tar.gz`, and `metaConfiguration.js` and `newColumnConfig.js` of `S_pneumoniae616`. The file `S_pneumoniae616/geneCluster/GC00000034_2_na_aln_reduced.fa` stays the archive's, because its S3 object holds records of other species

Together 32,464 GET and 66 LIST requests, about $0.014, and 2 GB of transfer, up to $0.18:

```bash
{
  printf '%s\n' index.json index.json.gz not_found.json not_found.json.gz
  jq -r '.Contents[]? | select(.LastModified | startswith("2022-09-09")) | .Key' "${listings}"/top-level/*.json
  printf '%s\n' \
    dataset/Bacillus_velezensis/geneCluster/d66a12ef3a081038_refined_na_aln.fa.gz \
    dataset/Citrobacter_portucalensis/geneCluster/9671524fd5b514c4_refined_aa_aln.fa.gz \
    dataset/Citrobacter_portucalensis/geneCluster/fba6d8c25065b05f_refined_tree.json.gz \
    dataset/Phaeobacter_inhibens/geneCluster/be2a1c40c0196838_refined_aa_aln_reduced.fa.gz \
    dataset/S_pneumoniae616/metaConfiguration.js dataset/S_pneumoniae616/newColumnConfig.js
  for dataset in S_pneumoniae616 P_marinus_meta Enterobacteriales Pseudomonadales Vibrionales; do
    printf "dataset/${dataset}/%s\n" all_gene_alignments.tar.gz core_gene_alignments.tar.gz
  done
} > s3-keys.txt
PANX_OUT="${objects}" dev/docker/run python3 -m panx_data fetch-s3 --keys s3-keys.txt --dest "${objects}"
aws s3 sync "s3://${bucket}/dataset/Neisseria_gonorrhoeae_ECDC/" "${objects}/dataset/Neisseria_gonorrhoeae_ECDC/"
```

### 4. Write the checksums and protect `src/`

```bash
PANX_OUT="${snapshot}/src" dev/docker/run bash -c \
  "cd '${snapshot}/src' && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum > MD5SUMS.tmp && mv MD5SUMS.tmp MD5SUMS"
chmod -R a-w "${snapshot}/src"
```

A release recipe names the SHA-256 of `src/MD5SUMS`, so a release builds only from the inputs that it names.

### 5. Check the inputs against the old bucket

`rebuild` writes the old layout from the inputs, and `verify` compares it with the listings. A release does not need this check again; it shows that the inputs hold the data of the old bucket. No AWS requests. The rebuild of all 270 datasets takes about 10 minutes with 8 jobs and writes 24,819,403 files (127.7 GB):

```bash
dev/docker/run just rebuild "${snapshot}" --jobs 16
dev/docker/run just verify "${snapshot}"
```

`verify` must report only known classes. For the bucket as listed on 2026-11-11, expect these counts:

| Class                                  |   Count | Meaning                                                                                   |
| -------------------------------------- | ------: | ----------------------------------------------------------------------------------------- |
| `missing-in-s3`                        |  16,054 | the upload left the object out; `out/` has the correct file                               |
| `s3-empty`, `s3-shorter`               |  5 + 12 | the upload broke the object; `out/` has the correct file                                  |
| `only-in-s3:stale-archive`             |     424 | an old `*_gene_alignments.tar.gz` that nothing links                                      |
| `differs:gzip-header-of-older-dataset` | 351,664 | a `.gz` copy of an older dataset that S3 holds with another gzip header; the same content |

A `differs` row means that S3 holds other bytes for a key. Copy that key into `src/s3/` as in step 3 (one GET request per key, $0.43 per 1M, plus its size in transfer), write `MD5SUMS` again, delete `out/pangenome.org-data-master/`, and rebuild.

### 6. Build and check a release

A recipe `releases/<release id>.toml` names the release id (its UTC creation time, such as `20261112T080000Z`), its title, the SHA-256 of `src/MD5SUMS`, and its datasets in display order. `release-build` writes the release folder into the data root `out/data.pangenome.org/` and lists it in its `releases.json`; `release-check` compares every release of the data root with the inputs and prints its findings. No AWS requests. The release of all 270 datasets of 2026-11-11 holds 2,162 objects and 11.38 GB:

```bash
dev/docker/run just release-build releases/<release id>.toml "${snapshot}" --jobs 8
dev/docker/run just release-check "${snapshot}/out/data.pangenome.org" "${snapshot}" | tee "${snapshot}/out/release-check.txt"
```

`release-build` refuses to start unless Python's zlib is the zlib of the pinned Python, because the deflate bytes depend on the zlib build. `release-check` exits 1 on a failing finding. `data-defects.md` lists the defects of the input data that the build handles.

### 7. Check the release and run the app on it

In the app repository `neherlab/pan-genome-visualization`, check the data root against the data format and serve the app with it (the app at port 8000, the data root at port 8001). No AWS requests:

```bash
PANX_DATA="${snapshot}/out/data.pangenome.org" dev/docker/run just verify-data
PANX_DATA="${snapshot}/out/data.pangenome.org" DEV_PUBLISH=8000,8001 dev/docker/run just serve
```

### 8. Publish the release

`release-publish` publishes a release of a local data root to the data host, the folder `<prefix>` of the data bucket, with the AWS CLI, one step at a time. Run each step first with `--dry-run`, which prints the requests, their count, and their cost, and sends nothing:

```bash
publish() { dev/docker/run just release-publish "$1" "${snapshot}/out/data.pangenome.org" <release id> --bucket <data bucket> --prefix data.pangenome.org "${@:2}"; }
publish upload --dry-run
```

1. `upload`: every object of the release folder with the headers of its manifest entry, then `manifest.json`. Each write sends the SHA-256 of the manifest, so S3 rejects other bytes, and is conditional on a missing key, so no object is replaced. A key that exists with the manifest's SHA-256 counts as done, so an interrupted upload resumes; a key with other bytes stops the upload, because a path of a release is written once. 2,162 PUT requests for the release of all 270 datasets, $0.012, plus one HEAD request per existing key; the release stores 11.38 GB, $0.28 per month
2. check: `just verify-data` of step 7 on the uploaded data root, then on the data host with `PANX_DATA=https://data.pangenome.org` and `--release <release id>`, about 2,700 HEAD and GET requests, below $0.01
3. `list`: add the release to `releases.json` of the data host and keep the default release; the release then opens with `?release=<release id>`. 1 HEAD, 1 GET, and 1 PUT request, below $0.0001
4. `promote`: make the release the default. Every visitor gets it within 10 minutes, because browsers and the CDN keep `releases.json` for at most 5 minutes each. 1 HEAD, 1 GET, and 1 PUT request, below $0.0001

On a data host without `releases.json`, `list` and `promote` create it with the release as its only entry and default. Every later write of `releases.json` is conditional on the ETag that the step read, so a second publish at the same time fails instead of dropping a release. A rollback is a `promote` of an earlier release. `unlist` takes a release that is not the default out of `releases.json`, the first step of its removal; keep its folder for at least 30 days after that, so open pages keep working ($0.28 for the 30 days).

`release-publish` refuses a release whose `manifest.json` records uncommitted changes of the tool, or another Python or zlib than `.config/mise.toml` pins, and a release folder whose files do not have the sizes of the manifest.

## Commands

```bash
python3 -m panx_data rebuild <snapshot>           # build <snapshot>/out/; --jobs N, --dataset NAME, --out DIR
python3 -m panx_data verify <snapshot>            # write out/deviations.tsv; exit 1 on an unknown difference
python3 -m panx_data fetch-s3 --keys FILE --dest <snapshot>/src/s3/pangenome.org-data-master
python3 -m panx_data release-build <recipe> <snapshot>       # --out DATA_ROOT, --jobs N, --dataset ID
python3 -m panx_data release-check <data root> <snapshot>    # exit 1 on a failing finding
python3 -m panx_data release-publish <step> <data root> <release id> --bucket B --prefix P [--dry-run]
```

`just --list` lists the recipes that run these commands, the tests, and the lints. `rebuild` resumes an interrupted run and refuses to overwrite a complete one; `release-build` refuses to overwrite a release folder. `fetch-s3` and `release-publish` use the AWS CLI with the caller's credentials; `fetch-s3` costs one GET request per object ($0.43 per 1M) plus the transfer of the objects.

## How the old bucket was filled

`panx_data/rebuild.py` describes the upload of 2022 that the rebuild repeats: the decompression and `.gz` copy of every file, the alignment zips, the August rule that kept the pipeline's original `.gz` file when the new copy had the same size, and the objects that S3 changed later.

`verify` sorts every difference into a class. The known classes describe the bucket, not the rebuild:

- `missing-in-s3`: the upload left the object out
- `s3-empty`, `s3-shorter`: the upload broke the object
- `only-in-s3:stale-archive`: an old `*_gene_alignments.tar.gz` that nothing links
- `only-in-s3:dataset-not-rebuilt`: a dataset that the archives lack
- `differs:gzip-header-of-older-dataset`: the 2022 upload compressed a file of an older dataset again, with another gzip header time

Any other class means that the rebuild does not reproduce the bucket.

## Tests and lints

```bash
dev/docker/run just check
```
