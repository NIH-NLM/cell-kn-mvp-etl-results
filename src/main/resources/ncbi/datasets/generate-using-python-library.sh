#!/usr/bin/env bash
PROJECT_NAME="ncbi-datasets-openapi"
openapi-generator-cli \
    generate \
    -g python \
    -i openapi3.docs.yaml \
    -o $PROJECT_NAME \
    --package-name "ncbi.datasets.openapi" \
    --additional-properties=pythonAttrNoneIfUnset=true,projectName=$PROJECT_NAME
