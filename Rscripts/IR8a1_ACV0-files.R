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
#path = "~/Research/projects/vikas-fly-1/IR8a1_ACV0/orig_data/"
path = "/ytmp/sgerber/projects/vikas/fly-data-1/fly-moves/IR8a1_ACV0/orig_data/"


xyFiles = c(
"130709video34_xypts_transformed.csv",
"130709video4_xypts_transformed.csv",
"130710video4_xypts_transformed.csv",
"130712video34_xypts_transformed.csv",
"130712video38_xypts_transformed.csv",
"130716video34_xypts_transformed.csv",
"130716video4_xypts_transformed.csv",
"130717video37_xypts_transformed.csv",
"130717video5_xypts_transformed.csv",
"130718video34_xypts_transformed.csv",
"130718video8_xypts_transformed.csv",
"130724video4_xypts_transformed.csv",
"130725video34_xypts_transformed.csv",
"130725video7_xypts_transformed.csv",
"130730video34_xypts_transformed.csv",
"130730video3_xypts_transformed.csv",
"130808video38_xypts_transformed.csv",
"130808video4_xypts_transformed.csv",
"130809video34_xypts_transformed.csv",
"130904video34_xypts_transformed.csv",
"130904video4_xypts_transformed.csv",
"130906video34_xypts_transformed.csv",
"130912video34_xypts_transformed.csv",
"130912video4_xypts_transformed.csv",
"130917video34_xypts_transformed.csv",
"130925video8_xypts_transformed.csv",
"131002video4_xypts_transformed.csv",
"131004video34_xypts_transformed.csv",
"131009video34_xypts_transformed.csv",
"131009video4_xypts_transformed.csv",
"131010video33_xypts_transformed.csv",
"131010video4_xypts_transformed.csv"
)




innerRimFiles <- c(
"130709_3transformedrimpoints-inner.csv",
"130709transformedrimpoints-inner.csv",
"130710transformedrimpoints-inner.csv",
"130712_3transformedrimpoints-inner.csv",
"130712_3transformedrimpoints-inner.csv",
"130716_3transformedrimpoints-inner.csv",
"130716transformedrimpoints-inner.csv",
"130717_3transformedrimpoints-inner.csv",
"130717transformedrimpoints-inner.csv",
"130718_3transformedrimpoints-inner.csv",
"130718transformedrimpoints-inner.csv",
"130724transformedrimpoints-inner.csv",
"130725_3transformedrimpoints-inner.csv",
"130725transformedrimpoints-inner.csv",
"130730transformedrimpoints-inner.csv",
"130730_3transformedrimpoints-inner.csv",
"130808_3transformedrimpoints-inner.csv",
"130808transformedrimpoints-inner.csv",
"130809_3transformedrimpoints-inner.csv",
"130904_3transformedrimpoints-inner.csv",
"130904transformedrimpoints-inner.csv",
"130906_3transformedrimpoints-inner.csv",
"130912_3transformedrimpoints-inner.csv",
"130912transformedrimpoints-inner.csv",
"130917_3transformedrimpoints-inner.csv",
"130925transformedrimpoints-inner.csv",
"131002transformedrimpoints-inner.csv",
"131004_3transformedrimpoints-inner.csv",
"131009_3transformedrimpoints-inner.csv",
"131009transformedrimpoints-inner.csv",
"131010_3transformedrimpoints-inner.csv",
"131010transformedrimpoints-inner.csv"
)



outerRimFiles <- c(
"130709_3transformedrimpoints-outer.csv",
"130709transformedrimpoints-outer.csv",
"130710transformedrimpoints-outer.csv",
"130712_3transformedrimpoints-outer.csv",
"130712_3transformedrimpoints-outer.csv",
"130716_3transformedrimpoints-outer.csv",
"130716transformedrimpoints-outer.csv",
"130717transformedrimpoints-outer.csv",
"130717_3transformedrimpoints-outer.csv",
"130718_3transformedrimpoints-outer.csv",
"130718transformedrimpoints-outer.csv",
"130724transformedrimpoints-outer.csv",
"130725_3transformedrimpoints-outer.csv",
"130725transformedrimpoints-outer.csv",
"130730transformedrimpoints-outer.csv",
"130730_3transformedrimpoints-outer.csv",
"130808_3transformedrimpoints-outer.csv",
"130808transformedrimpoints-outer.csv",
"130809_3transformedrimpoints-outer.csv",
"130904_3transformedrimpoints-outer.csv",
"130904transformedrimpoints-outer.csv",
"130906_3transformedrimpoints-outer.csv",
"130912_3transformedrimpoints-outer.csv",
"130912transformedrimpoints-outer.csv",
"130917_3transformedrimpoints-outer.csv",
"130925transformedrimpoints-outer.csv",
"131002transformedrimpoints-outer.csv",
"131004_3transformedrimpoints-outer.csv",
"131009_3transformedrimpoints-outer.csv",
"131009transformedrimpoints-outer.csv",
"131010_3transformedrimpoints-outer.csv",
"131010transformedrimpoints-outer.csv"
)

xyFiles =sprintf("%s%s", path, xyFiles)
innerRimFiles =sprintf("%s%s", path, innerRimFiles)
outerRimFiles =sprintf("%s%s", path, outerRimFiles)
