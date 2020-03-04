
##############
# From Samuel's code
## retry failed jobs
memory_increment_ratio  = float(config["memory_increment_ratio"]) if "memory_increment_ratio" in config.keys() else 1 
runtime_increment_ratio = float(config["runtime_increment_ratio"]) if "runtime_increment_ratio" in config.keys() else 1


## get incremental memory allocation when jobs fail
## 'startStr' is the first memory allocation in GB
## input is in GB; output in MB; default is 2GB, but can be changed by a rule

## get incremental memory allocation when jobs fail
## 'start' is the first memory allocation in GB (default 4GB)
## input is in GB; output is in MB;
## global variable memory_increment_ratio defines by how much (ratio) the memory is increased if not defined specifically
def get_memory_alloc(startStr, attempt, default=2):
	increment='{}_increment'.format(startStr)
	mem_start = int(config[startStr] if startStr in config.keys() else default)
	mem_incre = int(config[increment] if increment in config.keys() else memory_increment_ratio*mem_start)
	return int(1024 * ((attempt-1) * mem_incre + mem_start))
	
	
def convert_time(seconds): 
    day = seconds // (24 * 3600) 
    seconds = seconds % (24 * 3600) 
    hour = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60  
    return "%d-%02d:%02d:%02d" % (day, hour, minutes, seconds) 
	
	
## get incremental time allocation when jobs fail
## 'start' is the first time allocation in hours (default 12h)
## input is in hours; output is in minutes;
## global variable runtime_increment_ratio defines by how much (ratio) the time is increased if not defined specifically
def get_runtime_alloc(startStr, attempt, default=12):
	increment='{}_increment'.format(startStr)
	time_start = int(config[startStr] if startStr in config.keys() else default)
	time_incre = int(config[increment] if increment in config.keys() else runtime_increment_ratio*time_start)
	return int(60 * ((attempt-1) * time_incre + time_start))
	#return convert_time(60*60 * ((attempt-1) * time_incre + time_start))
	
	
## get the number of threads of the given parameter
def get_threads(param, default=1):
	return int(config[param] if param in config.keys() else default)
