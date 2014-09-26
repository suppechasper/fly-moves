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
path = "~/Research/projects/vikas-fly-1/Orco_ACV0/orig_data/"
#path = "/ytmp/sgerber/projects/vikas/fly-data-1/fly-moves/Orco_ACV0/orig_data/"


xyFiles = c(
"121011video10_xypts_transformed.csv",
"121012video5_xypts_transformed.csv",
"121017video3_xypts_transformed.csv",
"121018video3_xypts_transformed.csv",
"121019video3_xypts_transformed.csv",
"121019video7_xypts_transformed.csv",
"121023video2_xypts_transformed.csv",
"121024video4_xypts_transformed.csv",
"121024video8_xypts_transformed.csv",
"121025video4_xypts_transformed.csv",
"121025video6_xypts_transformed.csv",
"121026video2_xypts_transformed.csv",
"121026video5_xypts_transformed.csv",
"121030video3_xypts_transformed.csv",
"121030video6_xypts_transformed.csv",
"121031video4_xypts_transformed.csv",
"121031video7_xypts_transformed.csv",
"121101video3_xypts_transformed.csv",
"121106video3_xypts_transformed.csv",
"121106video6_xypts_transformed.csv",
"130226video7_xypts_transformed.csv",
"130304video36_xypts_transformed.csv",
"140326video33_xypts_transformed.csv",
"140401video35_xypts_transformed.csv",
"140402video35_xypts_transformed.csv",
"140408video33_xypts_transformed.csv",
"140409video32_xypts_transformed.csv",
"140411video33_xypts_transformed.csv",
"140415video2_xypts_transformed.csv",
"140415video33_xypts_transformed.csv"
)




innerRimFiles <- c(
"121011transformedrimpoints-inner.csv",
"121012transformedrimpoints-inner.csv",
"121017transformedrimpoints-inner.csv",
"121018transformedrimpoints-inner.csv",
"121019transformedrimpoints-inner.csv",
"121019transformedrimpoints-inner.csv",
"121023transformedrimpoints-inner.csv",
"121024transformedrimpoints-inner.csv",
"121024transformedrimpoints-inner.csv",
"121025transformedrimpoints-inner.csv",
"121025transformedrimpoints-inner.csv",
"121026transformedrimpoints-inner.csv",
"121026transformedrimpoints-inner.csv",
"121030transformedrimpoints-inner.csv",
"121030transformedrimpoints-inner.csv",
"121031transformedrimpoints-inner.csv",
"121031transformedrimpoints-inner.csv",
"121101transformedrimpoints-inner.csv",
"121106transformedrimpoints-inner.csv",
"121106transformedrimpoints-inner.csv",
"130226transformedrimpoints-inner.csv",
"130304transformedrimpoints-inner.csv",
"140326_3transformedrimpoints-inner.csv",
"140401_3transformedrimpoints-inner.csv",
"140402_3transformedrimpoints-inner.csv",
"140408_3transformedrimpoints-inner.csv",
"140409_3transformedrimpoints-inner.csv",
"140411_3transformedrimpoints-inner.csv",
"140415transformedrimpoints-inner.csv",
"140415_3transformedrimpoints-inner.csv"
)



outerRimFiles <- c(
"121011transformedrimpoints-outer.csv",
"121012transformedrimpoints-outer.csv",
"121017transformedrimpoints-outer.csv",
"121018transformedrimpoints-outer.csv",
"121019transformedrimpoints-outer.csv",
"121019transformedrimpoints-outer.csv",
"121023transformedrimpoints-outer.csv",
"121024transformedrimpoints-outer.csv",
"121024transformedrimpoints-outer.csv",
"121025transformedrimpoints-outer.csv",
"121025transformedrimpoints-outer.csv",
"121026transformedrimpoints-outer.csv",
"121026transformedrimpoints-outer.csv",
"121030transformedrimpoints-outer.csv",
"121030transformedrimpoints-outer.csv",
"121031transformedrimpoints-outer.csv",
"121031transformedrimpoints-outer.csv",
"121101transformedrimpoints-outer.csv",
"121106transformedrimpoints-outer.csv",
"121106transformedrimpoints-outer.csv",
"130226transformedrimpoints-outer.csv",
"130304transformedrimpoints-outer.csv",
"140326_3transformedrimpoints-outer.csv",
"140401_3transformedrimpoints-outer.csv",
"140402_3transformedrimpoints-outer.csv",
"140408_3transformedrimpoints-outer.csv",
"140409_3transformedrimpoints-outer.csv",
"140411_3transformedrimpoints-outer.csv",
"140415transformedrimpoints-outer.csv",
"140415_3transformedrimpoints-outer.csv"
)

xyFiles =sprintf("%s%s", path, xyFiles)
innerRimFiles =sprintf("%s%s", path, innerRimFiles)
outerRimFiles =sprintf("%s%s", path, outerRimFiles)

