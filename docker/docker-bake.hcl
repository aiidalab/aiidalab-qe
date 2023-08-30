# docker-bake.hcl for building QeApp images
group "default" {
  targets = ["qe"]
}

# The variable for the organization name, which usually is the same as the
# Docker Hub username or your GitHub organization name.
variable "ORGANIZATION" {
  default = "aiidalab"
}

target "qe" {
  tags = ["${ORGANIZATION}/qe:newly-baked"]
  context = "."
  contexts = {
    src = ".."
  }
  platforms = ["linux/amd64", "linux/arm64"]
}
