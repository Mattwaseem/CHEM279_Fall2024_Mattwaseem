FROM ubuntu:20.04

# Avoid interactive prompts during installation.
ENV DEBIAN_FRONTEND=noninteractive

# Install all necessary packages.
RUN apt-get update -y && \
    apt-get install -y \
    build-essential \
    cmake \
    nodejs \
    npm \
    libwebsocketpp-dev \
    curl \
    libboost-all-dev \
    liblapack-dev \
    libarmadillo-dev && \
    rm -rf /var/lib/apt/lists/*


# Step 5: Expose the port on which your Node.js app will run
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.3/install.sh | bash

# Copy the entrypoint file into the Docker image
COPY entrypoint.sh /entrypoint.sh

# Make the entrypoint script executable
RUN chmod +x /entrypoint.sh

# Define the entrypoint script that should be run
ENTRYPOINT ["/entrypoint.sh"]

EXPOSE 3000