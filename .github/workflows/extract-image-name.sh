# Extract image names together with their sha256 digests
# from the docker/bake-action metadata output.
# These together uniquely identify newly built images.

# The input to this script is a JSON string passed via BAKE_METADATA env variable
# Here's example input (trimmed to relevant bits):
# BAKE_METADATA: {
#    "base": {
#      "containerimage.descriptor": {
#        "mediaType": "application/vnd.docker.distribution.manifest.v2+json",
#        "digest": "sha256:8e57a52b924b67567314b8ed3c968859cad99ea13521e60bbef40457e16f391d",
#        "size": 6170,
#      },
#      "containerimage.digest": "sha256:8e57a52b924b67567314b8ed3c968859cad99ea13521e60bbef40457e16f391d",
#      "image.name": "ghcr.io/aiidalab/qe"
#    }
#  }
#
# Example output (real output is on one line):
#
# image="ghcr.io/aiidalab/qe@sha256:79a0f984b9e03b733304fda809ad3e8eec8416992ff334052d75da00cadb8f12"
# }
#
# This json output is later turned to environment variables using fromJson() GHA builtin
# (e.g. BUILD_MACHINE_IMAGE=ghcr.io/aiidalab/qe@sha256:8e57a52b...)
# and these are in turn read in the docker-compose.<target>.yml files for tests.

if [[ -z ${BAKE_METADATA-} ]];then
    echo "ERROR: Environment variable BAKE_METADATA is not set!"
    exit 1
fi

image=$(echo "${BAKE_METADATA}" | jq -c '. as $base | to_entries[] | [(.value."image.name"|split(",")[0]),(.value."containerimage.digest")]|join("@")' | tr -d '"')
echo "image=$image"
