#' azpro dataset
#'
#' Data come from the 1991 Arizona cardiovascular patient files. A subset of
#' the fields was selected to model the differential length of stay for
#' patients entering the hospital to receive one of two standard cardiovascular
#' procedures: CABG and PTCA. CABG is the standard acronym for Coronary
#' Artery Bypass Graft, where the flow of blood in a diseased or blocked
#' coronary artery or vein has been grafted to bypass the diseased sections.
#' PTCA, or Percutaneous Transluminal Coronary Angioplasty, is a method of
#' placing a balloon in a blocked coronary artery to open it to blood flow.
#' It is a much less severe method of treatment for those having coronary
#' blockage, with a corresponding reduction in risk.
#'
#' @format ## `azpro`
#' A data frame with 3589 observations on the following 6 variables:
#' \describe{
#'   \item{los}{length of hospital stay}
#'   \item{procedure}{1=CABG; 0=PTCA}
#'   \item{sex}{1=Male; 0=female}
#'   \item{admit}{1=Urgent/Emerg; 0=elective (type of admission)}
#'   \item{age75}{1=Age>75; 0=Age<=75}
#'   \item{hospital}{encrypted facility code (string)}
#' }
#' @source This dataset comes from the COUNT R package.
#' @author COUNT R package.
"azpro"
