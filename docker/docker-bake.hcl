# docker-bake.hcl for building QeApp images
group "default" {
  targets = ["qe"]
}

variable "QE_VERSION" {
}

variable "BASE_IMAGE" {
  default = "aiidalab/full-stack:latest"
}

variable "ORGANIZATION" {
  default = "aiidalab"
}

target "qe" {
  tags = ["${ORGANIZATION}/qe:newly-baked"]
  context = "."
  contexts = {
    src = ".."
    base-image = "docker-image://${BASE_IMAGE}"
  }
  args = {
    "QE_VERSION" = "${QE_VERSION}"
    "BASE_IMAGE" = "${BASE_IMAGE}"
  }
}
