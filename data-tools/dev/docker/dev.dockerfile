# syntax=docker/dockerfile:1
# check=experimental=all
#
# Development image of data-tools: the tools of .config/mise.lock (Python with
# its own zlib, uv, ruff, aws-cli, and the lints), and pigz 2.8 with zlib 1.3.1
# of Debian, which `rebuild` needs for byte-exact .gz copies of the 2022 upload.
# dev/docker/run builds it and runs every command as the host user, with
# HOME=/tmp/home.
FROM debian:13-slim@sha256:f6b08b09a05725f9431502d3b5714f5c4576e82c1a84fa1debe0c24a7e0c3ec4
SHELL ["bash", "-euxo", "pipefail", "-c"]

# The Debian packages come from the snapshot archive of the date of the base
# image (trixie-20261027-slim), so every build installs the same versions. The
# sources name only the snapshot, because apt would otherwise also read the
# current archive and take a newer version from it. Its Release files of
# trixie-updates and trixie-security have expired, so apt skips that check; it
# still checks their signatures. The base image already holds zlib1g, and apt
# keeps that version, so the digest fixes it. The build fails unless dpkg
# reports the versions that `rebuild` was checked with.
ARG DEBIAN_SNAPSHOT="20261027T000000Z"
ARG PIGZ_VERSION="2.8-1"
ARG ZLIB_VERSION="1:1.3.dfsg+really1.3.1-1+b1"
RUN set -euxo pipefail >/dev/null \
&& printf '%s\n' \
  "Types: deb" \
  "URIs: http://snapshot.debian.org/archive/debian/${DEBIAN_SNAPSHOT}" \
  "Suites: trixie trixie-updates" \
  "Components: main" \
  "Check-Valid-Until: no" \
  "Signed-By: /usr/share/keyrings/debian-archive-keyring.pgp" \
  "" \
  "Types: deb" \
  "URIs: http://snapshot.debian.org/archive/debian-security/${DEBIAN_SNAPSHOT}" \
  "Suites: trixie-security" \
  "Components: main" \
  "Check-Valid-Until: no" \
  "Signed-By: /usr/share/keyrings/debian-archive-keyring.pgp" \
  >"/etc/apt/sources.list.d/debian.sources"

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq -o Acquire::Retries=5 \
&& apt-get install --no-install-recommends --yes -qq -o Acquire::Retries=5 \
  ca-certificates \
  curl \
  git \
  "pigz=${PIGZ_VERSION}" \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

RUN set -euxo pipefail >/dev/null \
&& for expected in "pigz=${PIGZ_VERSION}" "zlib1g=${ZLIB_VERSION}"; do \
  actual="$(dpkg-query --show --showformat='${Version}' "${expected%%=*}")"; \
  if [[ "${actual}" != "${expected#*=}" ]]; then \
    printf 'Expected %s %s, found %s\n' "${expected%%=*}" "${expected#*=}" "${actual}" >&2; \
    exit 1; \
  fi; \
done

# mise installs every tool of .config/mise.toml at the URL and checksum of
# .config/mise.lock and links their executables into /usr/local/bin. The image
# then runs the same Python, and so the same zlib, as a host that runs
# `mise install`.
ENV MISE_DATA_DIR="/opt/mise"
ENV MISE_CACHE_DIR="/tmp/mise/cache"
ENV MISE_STATE_DIR="/tmp/mise/state"
COPY dev/docker/files/fetch dev/docker/files/checksums dev/docker/files/install-mise dev/docker/files/mise-version /
RUN set -euxo pipefail >/dev/null \
&& /install-mise "/usr/local/bin" \
&& rm /fetch /checksums /install-mise /mise-version
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/.config/
RUN set -euxo pipefail >/dev/null \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install \
&& bin_paths="$(mise -C "/tmp/mise/project" bin-paths)" \
&& printf '%s\n' "${bin_paths}" | while IFS= read -r dir; do \
  find -L "${dir}" -mindepth 1 -maxdepth 1 -type f -executable -exec ln -sf -t "/usr/local/bin" {} + || exit 1; \
done \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" \
&& python3 -c 'import sys, zlib; print(sys.version, zlib.ZLIB_RUNTIME_VERSION)' \
&& just --version \
&& uv --version

# dev/docker/run runs as the host user with HOME=/tmp/home.
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"
