import lsst.ts.schedulerConfig.scheduler_driver
assert type(config)==lsst.ts.schedulerConfig.scheduler_driver.SchedulerDriver, \
    'config is of type %s.%s instead of lsst.ts.schedulerConfig.scheduler_driver.SchedulerDriver' % (
        type(config).__module__, type(config).__name__)

config.startup_type = 'HOT'
config.startup_database = '/home/opsim/run_local/output/sextans_3732.db'
