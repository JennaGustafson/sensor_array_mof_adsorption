import yaml
import os.path

from redis import Redis
from rq import Queue

job_queue = None
redis_conn = None
redis_config = None

filename = 'redis.yaml'

if os.path.exists(filename):
    with open(filename, 'r') as yaml_file:
        redis_config = yaml.load(yaml_file)

    redis_conn = Redis(**redis_config['redis'])
    job_queue = Queue(redis_config['queue'], connection=redis_conn)
