
library(shinycssloaders)

ws <- function(ui) withSpinner(ui, type = 5,
                               color.background = "white")

navigation <- function() {
	includeHTML('../www/layout/navigation.html')
}

beginPage <- function() {
  includeHTML('../www/layout/beginpage.html')
}

endPage <- function() {
  includeHTML('../www/layout/endpage.html')
}