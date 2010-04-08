.onLoad <- function(a,b)
{
  require(futile)
  require(zoo)
  require(quantmod)
  .init()
}

.init <- function()
{
  logger <<- getLogger('tawny')
  tawny.options <<- OptionsManager('tawny.options',
    defaults=list(use.plots=FALSE))
}


