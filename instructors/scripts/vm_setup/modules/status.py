import json
from os import path

'''
This module monitors status of server initialization
'''


class JobStatus(object):
    def __init__(self, config, proj_path, hostname=None):
        self.status_filename = config.STATUS_FILE
        self.proj_path = proj_path
        self.status_file = None
        self.hostname = hostname
        self.init_status()
        if not self.host_exists():
            self.add_host()

    def host_exists(self):
        # get global status
        status_info = self.get_global_status()
        return self.hostname in status_info['status']

    def get_status(self, step):
        '''
        This function reads a json status file
        '''
        status = ''
        status_info = {}
        # load status file
        status_info = self.get_global_status()
        # check for the status of a given host and step
        if step in status_info['status'][self.hostname]:
            status = status_info['status'][self.hostname][step]
        else:
            status = 'not started'
        return status

    def update_status(self, step, step_status):
        '''
        This function updates the contents of status file
        '''
        # get global status
        status_info = self.get_global_status()

        # update the status for a step for a given host
        status_info['status'][self.hostname][step] = step_status
        self.update_global_status(status_info)

    def get_global_status(self):
        global_status = {}
        with open(self.status_file, 'r') as sf:
            global_status = json.load(sf)
        return global_status

    def update_global_status(self, global_status):
        with open(self.status_file, 'w') as sf:
            json.dump(global_status, sf, indent=2)

    def add_host(self):
        # get global status
        status_info = self.get_global_status()
        # add new hostname
        status_info['status'][self.hostname] = {}
        # update status info
        self.update_global_status(status_info)

    def init_status(self):
        '''
        This function initializes the status file in the project directory
        '''
        # if status file does not exist, create one
        status_info = {}
        self.status_file = path.join(self.proj_path, self.status_filename)
        if not path.exists(self.status_file):
            status_info['status'] = {}
            self.update_global_status(status_info)
