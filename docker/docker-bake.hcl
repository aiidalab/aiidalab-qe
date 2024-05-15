# docker-bake.hcl

variable "QE_VERSION" {
}

variable "BASE_IMAGE" {
}

variable "ORGANIZATION" {
  default = "aiidalab"
}

group "default" {
  targets = ["qe"]
}

target "qe" {
  tags = ["${REGISTRY}/${ORGANIZATION}/qe"]
  context = "."
  contexts = {
    src = ".."
    base-image = "docker-image://${BASE_IMAGE}"
  }
  args = {
    "QE_VERSION" = "${QE_VERSION}"
  }
}
