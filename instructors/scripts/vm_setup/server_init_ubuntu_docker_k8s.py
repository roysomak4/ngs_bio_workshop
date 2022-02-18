import argparse
import json
import os
from sys import exit

from config.server_init_config import ServerInitConf as srv_init_conf
from config.server_init_config import KubeInitConf as kube_init_conf
from modules.server_setup import ServerInit


def setup_server():
    # get commandline arguments
    args = parse_arguments()

    # create project. Returns absolute path
    proj_path = create_project(args.project_name)

    # sudoer file location
    sudoers = os.path.abspath('sudo/sudoers')

    # ssh config file location
    ssh_conf = os.path.abspath('ssh/sshd_config')

    # docker daemon.json file
    docker_daemon = os.path.abspath('docker/daemon.json')

    # parse the remote server info
    servers = get_remote_servers(args.server_list)

    for server in servers:
        print(f'Initiating server setup for {server["hostname"]}...')
        print('-------------------------------')
        server_init = ServerInit(
            server['hostname'],
            server['public'],
            server['private'],
            server['role'],
            args.non_root_user,
            srv_init_conf,
            sudoers,
            ssh_conf,
            docker_daemon,
            proj_path)
        server_init.setup_server()
        print(f'Remote server {server["hostname"]} initialized successfully.')
        print('-------------------------------')
    print('All servers initialized successfully!')

def create_project(proj_name):
    # remove leading and trailing spaces from project name
    # replace spaces with '_' in project name
    proj_name = proj_name.strip().replace(' ', '_')
    proj_path = os.path.abspath(f'projects/{proj_name}')
    # check if project exists
    if os.path.isdir(proj_path):
        print('Existing project found. Starting from where you left off in this project...')
        return proj_path
    else:
        os.mkdir(proj_path)
        return proj_path


def get_remote_servers(server_file):
    with open(server_file, 'r') as f:
        servers = json.load(f)
    return servers


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('project_name',
                        type=str,
                        help='Name of the project. Required field.')
    parser.add_argument('server_list',
                        type=str,
                        help='File containing IP addresses of remote servers. Required field.')
    parser.add_argument('non_root_user',
                        type=str,
                        help='username to create user without root previlages. Required field.')
    return parser.parse_args()


def str2bool(arg):
    '''
    This module handles boolean argument input
    Source: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    '''
    if isinstance(arg, bool):
        return arg
    elif arg.lower() in ['yes', 'y', '1', 'true', 't']:
        return True
    elif arg.lower() in ['no', 'n', '0', 'false', 'f']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == '__main__':
    setup_server()
