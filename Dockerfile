FROM ubuntu
MAINTAINER test-user
RUN apt update
CMD ["echo", "Hello World"]
