tag ?= latest

all:
	docker buildx build \
			--platform linux/amd64 \
    	--push \
    	--tag 195996028523.dkr.ecr.eu-west-1.amazonaws.com/sention-eval:${tag} \
			.

local:
	docker build --output=type=docker --progress=plain --tag sention-eval:${tag} .