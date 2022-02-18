from os import environ


class ServerInitConf:
    TARGET_OS = 'ubuntu'
    DOCKER_DDL_URL = 'https://download.docker.com/linux/ubuntu/gpg'
    DOCKER_BASH_COMPLETION = 'https://raw.githubusercontent.com/docker/machine/v0.16.0/contrib/completion/bash/docker-machine.bash'
    SUDOERS = 'sudo/sudoers'
    SSH_CONFIG = 'ssh/sshd_config'
    STATUS_FILE = 'status.json'
    # user has to create the SSH_KEY and KNOWN_HOSTS environment variable in their client machine. This can be added using the export command in .bash_profile (mac) or .bashrc (linux)
    CLIENT_SSH_KEY = environ.get('SSH_KEY')
    CLIENT_KNOWN_HOSTS = environ.get('KNOWN_HOSTS')
    PORT = 22
    TIMEOUT = 10


class KubeInitConf:
    KUBE_REPO_KEY = 'https://packages.cloud.google.com/apt/doc/apt-key.gpg'
    POD_NETWORK_CIDR = '10.244.0.0/16'  # as recommended by flannel documentation
    FLANNEL_MANIFEST = 'https://raw.githubusercontent.com/coreos/flannel/master/Documentation/kube-flannel.yml'
    HOSTS_FILE = '/etc/hosts'
