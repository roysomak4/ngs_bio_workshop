import logging
from os import path


class Logger(object):
    def __init__(self, project):
        self.logger = logging.getLogger(__name__)
        self.loglevel = logging.INFO
        self.formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(message)s')
        self.file_handler = logging.FileHandler(path.join(project, 'server_init.log'))
        self.stream_handler = logging.StreamHandler()

    def init_logger(self):
        self.file_handler.setFormatter(self.formatter)
        self.stream_handler.setFormatter(self.formatter)
        self.logger.setLevel(self.loglevel)
        self.logger.addHandler(self.file_handler)
        self.logger.addHandler(self.stream_handler)

    def log_info(self, msg):
        self.logger.info(msg)

    def log_error(self, err):
        self.logger.error(err)

    def log_warning(self, warning):
        self.logger.warn(warning)
