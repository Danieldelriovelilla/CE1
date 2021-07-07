STK.V.9.0
WrittenBy    StandardObjectCatalog
BEGIN Facility
Name        Kiruna_Station_KIR-2_13m
	BEGIN CentroidPosition
		DisplayCoords         Geodetic
		EcfLatitude           67.858428
		EcfLongitude          20.96688
		EcfAltitude           385.8
		DisplayAltRef         Ellipsoid
		AzElMask              AzElMaskFile: Kiruna_Station_KIR-2_13m.aem
	END CentroidPosition
BEGIN Extensions
	BEGIN Graphics
		BEGIN Graphics
			ShowAzElAtRangeMask       On
			MinDisplayRange           0.0
			MaxDisplayRange           1000000.0
			NumAzElAtRangeMaskSteps   10
		END Graphics
	END Graphics
	BEGIN AccessConstraints
			LineOfSight     IncludeIntervals
			AzElMask        IncludeIntervals
	END AccessConstraints
	BEGIN Desc
		ShortText    24
Kiruna Station KIR-2 13m
		LongText    402
Name:           Kiruna Station KIR-2 13m
Country:        Sweden
Location:       Kiruna, Salmijarvi
Status:         Active
Type:           GroundStation
Notes:          

Sources:                       http://www.esa.int/esaMI/Operations/SEMERDSMTWE_0.html
               http://track.sfo.jaxa.jp/en/index.html
               http://spot4.cnes.fr/spot4_fr/stat2ghz.htm
Last updated:   2010-12-15
	END Desc
END Extensions
END Facility
