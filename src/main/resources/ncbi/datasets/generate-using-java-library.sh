#!/usr/bin/env bash
PROJECT_NAME="ncbi-datasets-openapi"
java -jar openapi-generator-cli-7.2.0.jar generate \
    -g python \
    -i openapi3.docs.yaml \
    -o $PROJECT_NAME \
    --package-name "ncbi.datasets.openapi" \
    --additional-properties=pythonAttrNoneIfUnset=true,projectName=$PROJECT_NAME
