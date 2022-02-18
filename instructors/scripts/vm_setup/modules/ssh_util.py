import paramiko
import sh
import socket


class SSH_util(object):
    def __init__(self, hostname, pub_ip, priv_ip, app_config):
        # user provided arguments
        self.hostname = hostname
        self.pub_ip = pub_ip
        self.priv_ip = priv_ip
        # config options
        self.app_config = app_config
        self.client = None
        self.ssh_output = None
        self.ssh_err = None

    def add_ssh_key(self, user):
        print(f'Setting up ssh access to remote server {self.hostname}...')
        cmd = f'ssh-copy-id -i {self.app_config.CLIENT_SSH_KEY} -o StrictHostKeyChecking=no {user}@{self.pub_ip}'
        sh.bash('-c', cmd, _fg=True)

    def connect(self, user):
        connected = True
        # login to remote server and get a connection handle
        try:
            print(f'Negotiating connection to remote server {self.hostname}...')
            # init ssh client instance
            self.client = paramiko.SSHClient()
            # get private ssh key
            pkey = paramiko.RSAKey.from_private_key(open(self.app_config.CLIENT_SSH_KEY))
            # load known hosts list
            self.client.load_host_keys(self.app_config.CLIENT_KNOWN_HOSTS)
            # connect to remote server
            self.client.connect(hostname=self.pub_ip,
                                port=self.app_config.PORT,
                                username=user,
                                pkey=pkey,
                                timeout=self.app_config.TIMEOUT,
                                allow_agent=False,
                                look_for_keys=False)
            print(f'Connected to remote server {self.hostname}')
        except paramiko.AuthenticationException:
            print(f'Authentication failed for {self.hostname}, please verify credentials.')
            connected = False
        except paramiko.SSHException as se:
            print(f'Could not establish SSH connection with {self.hostname}: {se}')
            connected = False
        except socket.timeout as toe:
            print(f'Connection timeout to {self.hostname}: {toe}')
            connected = False
        except Exception as e:
            print(f'Exception occurred while connecting to {self.hostname}: {e}')
            connected = False
        return connected

    def exec_cmd(self, cmd):
        """Execute a command on the remote host. Return a tuple containing
        an integer status and a two strings, the first containing stdout
        and the second containing stderr from the command."""
        exec_result = True
        try:
            print(f'Executing {cmd}...')
            stdin, stdout, stderr = self.client.exec_command(cmd, timeout=self.app_config.TIMEOUT, get_pty=True)
            self.ssh_output = stdout.read().decode('utf-8')
            self.ssh_err = stderr.read().decode('utf-8')
            if self.ssh_err:
                print(f'Error executing command. {self.ssh_err}')
                exec_result = False
            else:
                print(self.ssh_output)
        except paramiko.SSHException as se:
            print(f'Failed to execute SSH command. {se}')
            exec_result = False
        except socket.timeout as toe:
            print(f'Command timed out. {toe}')
            exec_result = False
        finally:
            return exec_result

    def exec_commands(self, commands, user):
        """Execute multiple command on the remote host."""
        result = True
        if self.connect(user):
            for command in commands:
                result = self.exec_cmd(command)
                if not result:
                    break
        else:
            print('Could not establish SSH connection.')
            result = False
        self.client.close()
        return result

    def upload_file(self, user, localfile, remotefile):
        "This method uploads the file to remote server"
        upload_flag = True
        try:
            if self.connect(user):
                ftp_client = self.client.open_sftp()
                ftp_client.put(localfile, remotefile)
                ftp_client.close()
                self.client.close()
            else:
                print("Could not establish SSH connection")
                upload_flag = False
        except Exception as e:
            print(f'Unable to upload the file to the remote server {remotefile}')
            print(f'Error: {e}')
            upload_flag = False
        finally:
            ftp_client.close()
            self.client.close()
        return upload_flag
