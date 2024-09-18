#!MC 1410
$!OpenLayout  "./plotrms.lay"

$!XYLineAxis YDetail 1 {RangeMin = 0.0}
$!XYLineAxis YDetail 1 {RangeMax = 1.0}

$!AlterData 
  PartList = '"Adaptive-1.0-100" "Adaptive-2.0-100" "Adaptive-4.0-100" "Adaptive-8.0-100" "Adaptive-16.0-100"'
  IgnoreDivideByZero = Yes
  Equation = 'V1=1-0.1*V1'

$!LOOP 5
$!XYLineAxis XDetail 1 {RangeMin = 0.2}
$!XYLineAxis XDetail 1 {RangeMax = 0.45}
$!XYLineAxis XDetail 1 {RangeMin = 0.55}
$!XYLineAxis XDetail 1 {RangeMax = 0.8}
$!XYLineAxis XDetail 1 {Title{Text = 'Truncation correlation distance'}}
$!VarSet |I|=(|LOOP|)
$!LineMap [1-4]  Assign{Zone = |I|}
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ImageWidth = 2646
$!ExportSetup ExportFName = 'adaptive|I|.eps'
$!Export
  ExportRegion = AllFrames
$!ENDLOOP

$!LOOP 4
$!XYLineAxis XDetail 1 {RangeMin = 20}
$!XYLineAxis XDetail 1 {RangeMax = 160}
$!XYLineAxis XDetail 1 {Title{Text = 'Truncation Distance (grid cells)'}}
$!VarSet |I|=(|LOOP|+5)
$!LineMap [1-4]  Assign{Zone = |I|}
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = EPS
$!ExportSetup EPSPreviewImage{ImageType = None}
$!ExportSetup ImageWidth = 2646
$!ExportSetup ExportFName = 'distance|I|.eps'
$!Export
  ExportRegion = AllFrames
$!ENDLOOP
