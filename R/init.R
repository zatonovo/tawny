.onLoad <- function(a,b)
{
  require(futile)
  require(futile.logger)
  require(zoo)
  require(quantmod)
  .init()
}

.init <- function(loglevel=INFO)
{
  config_logger(threshold=loglevel)
  if (!exists('logger')) logger <<- getLogger('tawny')
  if (!exists('tawny.options'))
    tawny.options <<- OptionsManager('tawny.options',
      defaults=list(use.plots=FALSE))
}

