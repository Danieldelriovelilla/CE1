STK.V.9.0
WrittenBy    StandardObjectCatalog
BEGIN Facility
Name        Svalbard_STDN_S22S
	BEGIN CentroidPosition
		DisplayCoords         Geodetic
		EcfLatitude           78.232908
		EcfLongitude          15.381773
		EcfAltitude           481.103
		DisplayAltRef         Ellipsoid
		AzElMask              AzElMaskFile: Svalbard_STDN_S22S.aem
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
		ShortText    18
Svalbard STDN S22S
		LongText    399
Name:           Svalbard STDN S22S
Country:        Norway
Location:       Longyearbyen, Svalbard
Status:         Active
Type:           GroundStation
Notes:          NASA#1734 || S-band AzEl 7.3m || new, 12/2010

Sources:                       NASA Directory of Station Locations Mar 03 2011
Last updated:   2011-03-21Antennas:       
  Type    :ParabolicReflector
  Diameter: 7.3 [Meters]

	END Desc
END Extensions
END Facility
