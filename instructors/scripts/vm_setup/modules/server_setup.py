import sh
import sys

from .status import JobStatus
from .log_module import Logger


class ServerInit(object):
    def __init__(self, hostname, pub_ip, priv_ip, vm_role, user, app_config, sudoers, ssh_conf, docker_daemon, proj_path):
        # user provided arguments
        self.hostname = hostname
        self.pub_ip = pub_ip
        self.priv_ip = priv_ip
        self.vm_role = vm_role
        self.user = user
        # config options
        self.app_config = app_config
        self.server_conn = None
        self.sudoers = sudoers
        self.ssh_conf = ssh_conf
        self.daemon = docker_daemon
        self.status = JobStatus(app_config, proj_path, hostname)
        self.logger = Logger(proj_path)

    def setup_server(self):
        # initialize logging module
        self.logger.init_logger()
        # preseed server
        self.preseed_server()
        # configure sudo
        self.configure_sudo()
        self.logger = None

    def add_ssh_key(self, user):
        try:
            self.logger.log_info(f'Setting up ssh access to remote server {self.hostname} at {self.pub_ip}...')
            cmd = f'ssh-copy-id -i {self.app_config.CLIENT_SSH_KEY} -o StrictHostKeyChecking=no {user}@{self.pub_ip}'
            sh.bash('-c', cmd, _fg=True)
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            sys.exit(1)

    def preseed_server(self):
        try:
            # update status
            job_status = self.status.get_status('preseed')
            if job_status in ['not started', 'failed']:
                self.status.update_status('preseed', 'working')
                self.logger.log_info(
                    f'Preseeding {self.hostname} server, installing sudo if missing and creating non-root user {self.user}...')
                # create new user
                print('-------------------------------')
                self.logger.log_info('Adding new user bioseq and setting timezone...')
                cmd1 = (f'ssh -t root@{self.pub_ip} "'
                       f'useradd -m -U {self.user} -s /bin/bash'
                       f'echo -e \"bioseq123\nbioseq123\" | passwd {self.user}'
                       f'usermod -aG sudo {self.user}"')
                sh.bash("-c", cmd1, _fg=True)
                # update ssh config and add ssh key for newly created user
                print('-------------------------------')
                self.logger.log_info('updating ssh config...')
                cmd2 = (f'ssh -t root@{self.pub_ip} "'
                       'sed -i "s/^PasswordAuthentication no/PasswordAuthentication yes/" /etc/ssh/sshd_config'
                       'systemctl restart sshd'
                       f'cp -R /root/.ssh /home/{self.user}/'
                       f'chmod -R {self.user}:{self.user} /home/{self.user}/.ssh"')
                sh.bash("-c", cmd2, _fg=True)
                
                self.logger.log_info(f'{self.hostname} preseeded successfully')
                self.status.update_status('preseed', 'completed')
                print('-------------------------------')
            else:
                if job_status == 'completed':
                    self.logger.log_info(f'Preseeding completed for {self.hostname}. Skipping to sudo config...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('preseed', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('preseed', 'failed')
            sys.exit(1)

    def configure_sudo(self):
        try:
            job_status = self.status.get_status('sudo_config')
            if job_status in ['not started', 'failed']:
                self.status.update_status('sudo_config', 'working')
                self.logger.log_info(f'Configuring passwordless sudo for {self.hostname}...')
                # upload the sudoers files to tmp directory in remote server
                print('-------------------------------')
                self.logger.log_info('Uploading sudoers file...')
                cmd1 = f'scp {self.sudoers} root@{self.pub_ip}:/tmp/sudoers'
                sh.bash("-c", cmd1, _fg=True)
                # replace original sudoers with uploaded sudoers
                print('-------------------------------')
                self.logger.log_info('Updating sudoers on remote server...')
                cmd2 = (f'ssh -t root@{self.pub_ip} "'
                        'chmod 440 /tmp/sudoers && '
                        'chown root:root /tmp/sudoers && '
                        'cp /etc/sudoers /tmp/sudoers.ori && '
                        'mv /tmp/sudoers /etc/"')
                sh.bash("-c", cmd2, _fg=True)
                self.logger.log_info(f'Passwordless sudo configured for {self.hostname}.')
                self.logger.log_info('Original sudoers file backed up at /tmp/sudoers.ori')
                print('-------------------------------')
                self.status.update_status('sudo_config', 'completed')
            else:
                if job_status == 'completed':
                    self.logger.log_info(f'Sudo is already configured for {self.hostname}. Skipping to secure ssh...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('sudo_config', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('sudo_config', 'failed')
            sys.exit(1)

    def configure_secure_ssh(self):
        try:
            job_status = self.status.get_status('secure_ssh')
            if job_status in ['not started', 'failed']:
                self.status.update_status('secure_ssh', 'working')
                self.logger.log_info(f'Configuring secure SSH for {self.hostname}...')
                # upload sshd_config file to remote server
                print('-------------------------------')
                self.logger.log_info('Uploading sshd_config to remote server...')
                cmd1 = f'scp {self.ssh_conf} {self.user}@{self.pub_ip}:/tmp/sshd_config'
                sh.bash("-c", cmd1, _fg=True)
                # update sshd_config settings
                print('-------------------------------')
                self.logger.log_info('Updating sshd_config on remote server...')
                cmd2 = (f'ssh -t {self.user}@{self.pub_ip} "'
                        'sudo chown root:root /tmp/sshd_config && '
                        'sudo mv /tmp/sshd_config /etc/ssh/ && '
                        'sudo systemctl restart ssh"')
                sh.bash("-c", cmd2, _fg=True)
                self.logger.log_info(f'Secure SSH configured for {self.hostname}.')
                self.status.update_status('secure_ssh', 'completed')
                print('-------------------------------')
            else:
                if job_status == 'completed':
                    self.logger.log_info(
                        f'Secure SSH is already configured for {self.hostname}. Skipping to firewall config...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('secure_ssh', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('secure_ssh', 'failed')
            sys.exit(1)

    def configure_firewall(self):
        try:
            job_status = self.status.get_status('firewall_config')
            print(job_status)
            if job_status in ['not started', 'failed']:
                self.status.update_status('firewall_config', 'working')
                self.logger.log_info(f'Configuring firewall using ufw on {self.hostname}...')
                cmd = (f'ssh -t {self.user}@{self.pub_ip} "'
                       'sudo ufw default deny incoming && '
                       'sudo ufw default allow outgoing && '
                       'sudo ufw allow 22 && '
                       'sudo ufw allow 443/tcp && '
                       'sudo ufw allow 80/tcp && '
                       'sudo ufw allow 6443/tcp && '
                       'sudo ufw allow 2379:2380/tcp && '
                       'sudo ufw allow 10248:10259/tcp && '
                       'sudo ufw allow 8472/udp && '
                       'sudo ufw allow 8472/tcp && '
                       'sudo ufw allow 8285/udp && '
                       'sudo ufw allow 8285/tcp && '
                       'yes | sudo ufw enable && '
                       'sudo ufw reload && '
                       'sudo ufw status"')
                sh.bash("-c", cmd, _fg=True)
                if self.vm_role == 'worker':
                    cmd2 = (f'ssh -t {self.user}@{self.pub_ip} "'
                            'sudo ufw allow 30000:32767/tcp"'
                            )
                    sh.bash("-c", cmd2, _fg=True)
                self.logger.log_info(f'Firewall configured and activated on {self.hostname}.')
                self.status.update_status('firewall_config', 'completed')
                print('-------------------------------')
            else:
                if job_status == 'completed':
                    self.logger.log_info(
                        f'Firewall is already configured for {self.hostname}. Skipping to docker installation...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('firewall_config', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('firewall_config', 'failed')
            sys.exit(1)

    def install_docker(self):
        try:
            job_status = self.status.get_status('docker')
            if job_status in ['not started', 'failed']:
                self.status.update_status('docker', 'working')
                self.logger.log_info(f'Installing Docker on {self.hostname}...')
                cmd1 = (
                    f'ssh -t {self.user}@{self.pub_ip} "'
                    'sudo apt-get update && '
                    'sudo apt-get install -y -q --no-install-recommends libapparmor1 aufs-tools ca-certificates libltdl7 apt-transport-https curl gnupg2 software-properties-common nfs-common conntrack && '
                    f'curl -fsSL {self.app_config.DOCKER_DDL_URL} | sudo apt-key add -"')
                sh.bash("-c", cmd1, _fg=True)
                print('-------------------------------')
                cmd2 = f"ssh -t {self.user}@{self.pub_ip} 'sudo add-apt-repository \"deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable\"'"
                sh.bash("-c", cmd2, _fg=True)
                print('-------------------------------')
                cmd3 = (
                    f'ssh -t {self.user}@{self.pub_ip} "'
                    'sudo apt-get update && '
                    'sudo apt-get install -y -q docker-ce && '
                    f'sudo usermod -aG docker {self.user} && '
                    'sudo systemctl enable --now docker && '
                    'docker --version && '
                    f'sudo curl -L {self.app_config.DOCKER_BASH_COMPLETION} -o /etc/bash_completion.d/docker-machine"')
                sh.bash("-c", cmd3, _fg=True)
                self.logger.log_info('Uploading docker daemon.json...')
                cmd4 = f'scp {self.daemon} {self.user}@{self.pub_ip}:/tmp/daemon.json'
                sh.bash("-c", cmd4, _fg=True)
                # update docker config
                print('-------------------------------')
                self.logger.log_info('Configuring docker to used systemd...')
                cmd5 = (f'ssh -t {self.user}@{self.pub_ip} "'
                        'sudo chown root:root /tmp/daemon.json && '
                        'sudo mv /tmp/daemon.json /etc/docker/ && '
                        'sudo mkdir -p /etc/systemd/system/docker.service.d && '
                        'sudo systemctl daemon-reload && '
                        'sudo systemctl restart docker"')
                sh.bash("-c", cmd5, _fg=True)
                self.logger.log_info(f'Docker daemon configured to use systemd for {self.hostname}.')
                self.logger.log_info(f'Installed Docker on {self.hostname}.')
                self.status.update_status('docker', 'completed')
                print('-------------------------------')
            else:
                if job_status == 'completed':
                    self.logger.log_info(f'Docker installed for {self.hostname}. Skipping...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('docker', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('docker', 'failed')
            sys.exit(1)

def install_ngs_apps(self):
        try:
            job_status = self.status.get_status('ngs_apps')
            if job_status in ['not started', 'failed']:
                self.status.update_status('docker', 'working')
                self.logger.log_info(f'Installing Docker on {self.hostname}...')
                cmd1 = (
                    f'ssh -t {self.user}@{self.pub_ip} "'
                    'sudo apt-get update && '
                    'sudo apt-get install -y -q --no-install-recommends libapparmor1 aufs-tools ca-certificates libltdl7 apt-transport-https curl gnupg2 software-properties-common nfs-common conntrack && '
                    f'curl -fsSL {self.app_config.DOCKER_DDL_URL} | sudo apt-key add -"')
                sh.bash("-c", cmd1, _fg=True)
                print('-------------------------------')
                cmd2 = f"ssh -t {self.user}@{self.pub_ip} 'sudo add-apt-repository \"deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable\"'"
                sh.bash("-c", cmd2, _fg=True)
                print('-------------------------------')
                cmd3 = (
                    f'ssh -t {self.user}@{self.pub_ip} "'
                    'sudo apt-get update && '
                    'sudo apt-get install -y -q docker-ce && '
                    f'sudo usermod -aG docker {self.user} && '
                    'sudo systemctl enable --now docker && '
                    'docker --version && '
                    f'sudo curl -L {self.app_config.DOCKER_BASH_COMPLETION} -o /etc/bash_completion.d/docker-machine"')
                sh.bash("-c", cmd3, _fg=True)
                self.logger.log_info('Uploading docker daemon.json...')
                cmd4 = f'scp {self.daemon} {self.user}@{self.pub_ip}:/tmp/daemon.json'
                sh.bash("-c", cmd4, _fg=True)
                # update docker config
                print('-------------------------------')
                self.logger.log_info('Configuring docker to used systemd...')
                cmd5 = (f'ssh -t {self.user}@{self.pub_ip} "'
                        'sudo chown root:root /tmp/daemon.json && '
                        'sudo mv /tmp/daemon.json /etc/docker/ && '
                        'sudo mkdir -p /etc/systemd/system/docker.service.d && '
                        'sudo systemctl daemon-reload && '
                        'sudo systemctl restart docker"')
                sh.bash("-c", cmd5, _fg=True)
                self.logger.log_info(f'Docker daemon configured to use systemd for {self.hostname}.')
                self.logger.log_info(f'Installed Docker on {self.hostname}.')
                self.status.update_status('docker', 'completed')
                print('-------------------------------')
            else:
                if job_status == 'completed':
                    self.logger.log_info(f'Docker installed for {self.hostname}. Skipping...')
        except sh.ErrorReturnCode as e:
            self.logger.log_error(f'Error: {str(e.stderr)}')
            self.status.update_status('docker', 'failed')
            sys.exit(2)
        except:
            self.logger.log_error(str(sys.exc_info()))
            self.status.update_status('docker', 'failed')
            sys.exit(1)
