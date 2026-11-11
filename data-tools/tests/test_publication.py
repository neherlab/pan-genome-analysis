import base64
import hashlib

import pytest

from panx_data.publication import (
    Release,
    head_matches,
    refusals,
    releases_requests,
    summary,
    updated_releases,
    upload_requests,
)

ID = "20261112T080000Z"
OLDER = "20250101T000000Z"
SHA = hashlib.sha256(b"pack").hexdigest()
OBJECTS = [
    {
        "path": "Testus/clusters.json",
        "size": 10,
        "sha256": hashlib.sha256(b"clusters").hexdigest(),
        "contentType": "application/json",
        "contentEncoding": "gzip",
        "cacheControl": "public, max-age=31536000, immutable",
        "contentDisposition": None,
    },
    {
        "path": "Testus/gene_clusters.tar.gz",
        "size": 20,
        "sha256": SHA,
        "contentType": "application/gzip",
        "contentEncoding": None,
        "cacheControl": "public, max-age=31536000, immutable",
        "contentDisposition": 'attachment; filename="Testus_gene_clusters.tar.gz"',
    },
]
BUILD = {
    "repository": "https://github.com/neherlab/pan-genome-analysis",
    "commit": "a" * 40,
    "uncommitted_changes": False,
    "python": "3.15.2",
    "zlib": "1.3.2",
    "snapshot": "b" * 64,
    "recipe": "c" * 64,
}
SIZES = {"Testus/clusters.json": 10, "Testus/gene_clusters.tar.gz": 20}


def release(**build: object) -> Release:
    manifest = {"format": 1, "release": ID, "build": {**BUILD, **build}, "objects": OBJECTS}
    return Release(ID, {"id": ID, "format": 1, "title": "panX data"}, manifest, "d" * 64, 30)


def entry(release_id: str) -> dict[str, object]:
    return {"id": release_id, "format": 1, "title": "panX data"}


def test_upload_sends_every_object_with_its_manifest_headers_and_the_manifest_last():
    requests = upload_requests(release(), "data.pangenome.org")

    assert [request.key for request in requests] == [
        f"data.pangenome.org/{ID}/Testus/clusters.json",
        f"data.pangenome.org/{ID}/Testus/gene_clusters.tar.gz",
        f"data.pangenome.org/{ID}/manifest.json",
    ]
    pack = requests[1].argv("bucket")
    assert pack[:7] == ["aws", "s3api", "put-object", "--bucket", "bucket", "--key", requests[1].key]
    assert "--content-encoding" not in pack
    for option, value in [
        ("--content-type", "application/gzip"),
        ("--cache-control", "public, max-age=31536000, immutable"),
        ("--content-disposition", 'attachment; filename="Testus_gene_clusters.tar.gz"'),
        ("--checksum-algorithm", "SHA256"),
        ("--checksum-sha256", base64.b64encode(bytes.fromhex(SHA)).decode()),
        ("--if-none-match", "*"),
    ]:
        assert pack[pack.index(option) + 1] == value
    manifest = requests[2].argv("bucket")
    assert manifest[manifest.index("--content-encoding") + 1] == "gzip"
    assert manifest[manifest.index("--checksum-sha256") + 1] == base64.b64encode(bytes.fromhex("d" * 64)).decode()


def test_upload_key_has_no_double_slash_for_a_prefix_with_slashes():
    assert upload_requests(release(), "/data.pangenome.org/")[0].key.startswith(f"data.pangenome.org/{ID}/")


def test_a_release_of_the_pinned_runtime_with_all_files_has_no_refusal():
    assert refusals(release(), SIZES, "3.15.2", "1.3.2") == []


@pytest.mark.parametrize(
    ("build", "reason"),
    [
        ({"uncommitted_changes": True}, "uncommitted"),
        ({"python": "3.13.0"}, "Python 3.13.0"),
        ({"zlib": "1.3.1"}, "zlib 1.3.1"),
    ],
)
def test_refuses_a_release_of_a_changed_tool_or_another_runtime(build, reason):
    [refusal] = refusals(release(**build), SIZES, "3.15.2", "1.3.2")
    assert reason in refusal


def test_refuses_a_release_with_a_missing_or_resized_file():
    found = refusals(release(), {"Testus/clusters.json": 11}, "3.15.2", "1.3.2")
    assert len(found) == 2
    assert "Testus/clusters.json" in found[0] and "11" in found[0]
    assert "Testus/gene_clusters.tar.gz" in found[1]


def test_refuses_a_manifest_of_another_release():
    manifest = {**release().manifest, "release": OLDER}
    other = Release(ID, entry(ID), manifest, "d" * 64, 30)
    assert any(OLDER in refusal for refusal in refusals(other, SIZES, "3.15.2", "1.3.2"))


def test_an_existing_key_counts_as_done_only_with_the_manifest_checksum():
    assert head_matches({"ChecksumSHA256": base64.b64encode(bytes.fromhex(SHA)).decode()}, SHA)
    assert not head_matches({"ChecksumSHA256": base64.b64encode(bytes.fromhex("e" * 64)).decode()}, SHA)
    assert not head_matches({}, SHA)


@pytest.mark.parametrize("step", ["list", "promote"])
def test_the_first_release_creates_releases_json_as_its_only_entry_and_default(step):
    assert updated_releases(None, step, entry(ID)) == {"format": 1, "default": ID, "releases": [entry(ID)]}


def test_list_adds_the_release_newest_first_and_keeps_the_default():
    current = {"format": 1, "default": OLDER, "releases": [entry(OLDER)]}
    assert updated_releases(current, "list", entry(ID)) == {
        "format": 1,
        "default": OLDER,
        "releases": [entry(ID), entry(OLDER)],
    }


def test_list_of_a_listed_release_writes_nothing():
    current = {"format": 1, "default": ID, "releases": [entry(ID)]}
    assert updated_releases(current, "list", entry(ID)) is None


def test_promote_then_promote_back_is_a_rollback():
    listed = {"format": 1, "default": OLDER, "releases": [entry(ID), entry(OLDER)]}
    promoted = updated_releases(listed, "promote", entry(ID))
    assert promoted == {**listed, "default": ID}
    assert updated_releases(promoted, "promote", entry(OLDER)) == listed


def test_promote_of_an_unlisted_release_lists_it():
    current = {"format": 1, "default": OLDER, "releases": [entry(OLDER)]}
    assert updated_releases(current, "promote", entry(ID)) == {
        "format": 1,
        "default": ID,
        "releases": [entry(ID), entry(OLDER)],
    }


def test_unlist_removes_a_release_that_is_not_the_default():
    current = {"format": 1, "default": OLDER, "releases": [entry(ID), entry(OLDER)]}
    assert updated_releases(current, "unlist", entry(ID)) == {
        "format": 1,
        "default": OLDER,
        "releases": [entry(OLDER)],
    }
    assert updated_releases({**current, "releases": [entry(OLDER)]}, "unlist", entry(ID)) is None


def test_unlist_refuses_the_default_and_a_data_host_without_releases_json():
    with pytest.raises(ValueError, match="default"):
        updated_releases({"format": 1, "default": ID, "releases": [entry(ID)]}, "unlist", entry(ID))
    with pytest.raises(ValueError, match="releases.json"):
        updated_releases(None, "unlist", entry(ID))


def test_releases_json_write_is_conditional_on_the_read():
    head, get, put = releases_requests("data.pangenome.org", ID, "promote")
    assert (head.operation, head.key) == ("head-object", f"data.pangenome.org/{ID}/manifest.json")
    assert (get.operation, get.key) == ("get-object", "data.pangenome.org/releases.json")
    argv = put.argv("bucket")
    assert put.operation == "put-object"
    assert argv[argv.index("--cache-control") + 1] == "public, max-age=300"
    assert "--if-match" in argv
    assert [request.operation for request in releases_requests("p", ID, "unlist")] == ["get-object", "put-object"]


def test_summary_counts_requests_and_their_cost():
    text = summary(upload_requests(release(), "p"))
    assert "3 requests" in text and "3 put-object" in text
    assert "$0.0000162" in text
