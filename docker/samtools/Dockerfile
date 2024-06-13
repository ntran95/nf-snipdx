#################################################################
# Dockerfile
#
# this Dockerfile creates an image from the biocontainters
# samtools docker image that can be used in GCP GLS
#################################################################

# Pull the appropriate samtools image
FROM biocontainers/samtools:v1.9-4-deb_cv1 as samtools

# using this as the new base image
FROM  ubuntu:16.04

# copy all data over
COPY --from=samtools / /

# add /usr/bin to PATH in case it isn't there
ENV PATH /usr/bin:$PATH

CMD ["/bin/bash"]
