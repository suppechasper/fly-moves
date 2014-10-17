#----- Extract data for analysis from raw data----


#---- Experimental setup information ----
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels


#Set of files from which data is extracted

#Set path to were files are located, if script is run in folder were data is
#located the current setting will work.
path = "~/Research/projects/vikas-fly-1/WT_ACV0/orig-data/"
#path = "/ytmp/sgerber/projects/vikas/fly-data-1/fly-moves/WT_ACV0/orig-data/"


xyFiles = c(
"120202video12_xypts_transformed.csv",
"120202video18_xypts_transformed.csv",
"120202video6_xypts_transformed.csv",
"120209video6_xypts_transformed.csv",
"120218video13_xypts_transformed.csv",
"120218video7_xypts_transformed.csv",
"120316video5_xypts_transformed.csv",
"120323video6_xypts_transformed.csv",
"120328video7_xypts_transformed.csv",
"120418video7_xypts_transformed.csv",
"120427video7_xypts_transformed.csv",
"120503video6_xypts_transformed.csv",
"120517video16_xypts_transformed.csv",
"120530video8_xypts_transformed.csv",
"120605video38_xypts_transformed.csv",
"120607video6_xypts_transformed.csv",
"120612video36_xypts_transformed.csv",
"120613video6_xypts_transformed.csv",
"120614video8_xypts_transformed.csv",
"120808video13_xypts_transformed.csv",
"120809video7_xypts_transformed.csv",
"120828video6_xypts_transformed.csv",
"120829video1_xypts_transformed.csv",
"120911video8_xypts_transformed.csv",
"120913video7_xypts_transformed.csv",
"120920video8_xypts_transformed.csv",
"120925video6_xypts_transformed.csv",
"120928video6_xypts_transformed.csv",
"130213video8_xypts_transformed.csv",
"130214video9_xypts_transformed.csv"
)




innerRimFiles <- c(
"120202transformedrimpoints.matinner_transformed.csv",
"120202transformedrimpoints.matinner_transformed.csv",
"120202transformedrimpoints.matinner_transformed.csv",
"120209transformedrimpoints.matinner_transformed.csv",
"120218transformedrimpoints.matinner_transformed.csv",
"120218transformedrimpoints.matinner_transformed.csv",
"120316transformedrimpoints.matinner_transformed.csv",
"120323transformedrimpoints.matinner_transformed.csv",
"120328transformedrimpoints.matinner_transformed.csv",
"120418transformedrimpoints.matinner_transformed.csv",
"120427transformedrimpoints.matinner_transformed.csv",
"120503transformedrimpoints.matinner_transformed.csv",
"120517transformedrimpoints.matinner_transformed.csv",
"120530transformedrimpoints.matinner_transformed.csv",
"120605transformedrimpoints.matinner_transformed.csv",
"120607transformedrimpoints.matinner_transformed.csv",
"120612transformedrimpoints.matinner_transformed.csv",
"120613transformedrimpoints.matinner_transformed.csv",
"120614transformedrimpoints.matinner_transformed.csv",
"120808transformedrimpoints.matinner_transformed.csv",
"120809transformedrimpoints.matinner_transformed.csv",
"120828transformedrimpoints.matinner_transformed.csv",
"120829transformedrimpoints.matinner_transformed.csv",
"120911transformedrimpoints.matinner_transformed.csv",
"120913transformedrimpoints.matinner_transformed.csv",
"120920transformedrimpoints.matinner_transformed.csv",
"120925transformedrimpoints.matinner_transformed.csv",
"120928transformedrimpoints.matinner_transformed.csv",
"130213transformedrimpoints.matinner_transformed.csv",
"130214transformedrimpoints.matinner_transformed.csv"
)



outerRimFiles <- c(
"120202transformedrimpoints.matouter_transformed.csv",
"120202transformedrimpoints.matouter_transformed.csv",
"120202transformedrimpoints.matouter_transformed.csv",
"120209transformedrimpoints.matouter_transformed.csv",
"120218transformedrimpoints.matouter_transformed.csv",
"120218transformedrimpoints.matouter_transformed.csv",
"120316transformedrimpoints.matouter_transformed.csv",
"120323transformedrimpoints.matouter_transformed.csv",
"120328transformedrimpoints.matouter_transformed.csv",
"120418transformedrimpoints.matouter_transformed.csv",
"120427transformedrimpoints.matouter_transformed.csv",
"120503transformedrimpoints.matouter_transformed.csv",
"120517transformedrimpoints.matouter_transformed.csv",
"120530transformedrimpoints.matouter_transformed.csv",
"120605transformedrimpoints.matouter_transformed.csv",
"120607transformedrimpoints.matouter_transformed.csv",
"120612transformedrimpoints.matouter_transformed.csv",
"120613transformedrimpoints.matouter_transformed.csv",
"120614transformedrimpoints.matouter_transformed.csv",
"120808transformedrimpoints.matouter_transformed.csv",
"120809transformedrimpoints.matouter_transformed.csv",
"120828transformedrimpoints.matouter_transformed.csv",
"120829transformedrimpoints.matouter_transformed.csv",
"120911transformedrimpoints.matouter_transformed.csv",
"120913transformedrimpoints.matouter_transformed.csv",
"120920transformedrimpoints.matouter_transformed.csv",
"120925transformedrimpoints.matouter_transformed.csv",
"120928transformedrimpoints.matouter_transformed.csv",
"130213transformedrimpoints.matouter_transformed.csv",
"130214transformedrimpoints.matouter_transformed.csv"
)

xyFiles =sprintf("%s%s", path, xyFiles)
innerRimFiles =sprintf("%s%s", path, innerRimFiles)
outerRimFiles =sprintf("%s%s", path, outerRimFiles)

