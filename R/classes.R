#

setClass("tsLib", representation(
	Name    = "character",
	RI      = "numeric",
	medRI   = "numeric",
	RIdev   = "matrix",
	selMass = "list",
	topMass = "list",
	libData = "data.frame",
	spectra = "list"
))

setClass("tsRim", representation(
	limits   = "matrix",
	standard = "numeric",
	mass     = "numeric")
)

setClass("tsSample", representation(
	Names    = "character",
	CDFfiles = "character",
	RIfiles  = "character",
	CDFpath  = "character",
	RIpath   = "character",
	days     = "character",
	data     = "data.frame")
)

setClass("tsMSdata", representation(
	RI        = "matrix",
	RT        = "matrix",
	Intensity = "matrix"
))

setClass("tsProfile", representation(
 	"tsMSdata",
 	info = "data.frame"
))

