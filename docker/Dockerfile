from debian:stable-slim
RUN apt update && apt install -y procps nano python3.9 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 1
COPY sentieon-genomics-202112.01 /sentieon
ENV SENTIEON_PYTHON=/usr/bin/python3