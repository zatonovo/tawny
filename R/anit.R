.onLoad <- function(a,b)
{
  require(futile)
  require(futile.logger)
  require(zoo)
  require(quantmod)
  .init()
}

.init <- function()
{
  if (!exists('logger')) logger <<- getLogger('tawny')
  if (!exists('tawny.options'))
    tawny.options <<- OptionsManager('tawny.options',
      defaults=list(use.plots=FALSE))
}

.init()
