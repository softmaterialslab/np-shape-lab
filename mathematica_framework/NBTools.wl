(* ::Package:: *)

(* Practical Functions:  workflow optimization, styling, etc. *)

(* Function to open this file, if you have it in the designated place: *)
toolFile:=SystemOpen@FileNameJoin[{$UserDocumentsDirectory,"Research Data\\NBTools.wl"}]

(* Function to open the current directory in your System default viewer (e.g. Windows Explorer): *)
curDir:=SystemOpen[Directory[]]

(* Function to open the directory containing the current Notebook: *)
noteDir:=SystemOpen[NotebookDirectory[]]

(* Function to set the directory to various important ones: *)
docDir:=Module[{},SetDirectory[$HomeDirectory];SetDirectory["Documents"];]
testDir:=Module[{},docDir; SetDirectory["Research Data\\Test_Bed"];]
scratchDir:=SetDirectory["E:\\Scratch"]

(* Function that sorts strings (typically files) by their numerical components (and deprecated version sorting by only first number): *)
	(*sortFiles[files_]:=files[[Ordering@PadRight@StringSplit[files,x:DigitCharacter..:>FromDigits@x]]]*)
sortFiles[files_]:=With[{numCount=Length@ToExpression[StringCases[files[[1]],NumberString]]},
SortBy[files,ToExpression[StringCases[#1,NumberString][[1;;numCount]]]&]]

(* Function to style a string with certain formatting: *)
myStyle[text_,othOpts_:{}]:=Style[text,FontFamily->"Times New Roman",22,Black,othOpts]

(* Label plots: *)
labelPlot[plot_, {xLabel_, yLabel_}]:=Labeled[plot, myStyle/@{xLabel,Rotate[yLabel,90\[Degree]]}, {Bottom, Left}]

(* Set the front end options to use black (not gray) text of my chosen font, at least in Graphics(may not apply to plot axes etc): *)
SetOptions[$FrontEndSession,GraphicsBoxOptions->{AxesStyle->Directive[FontFamily->"Times New Roman",Black]}]

simParams[directory_]:=Internal`StringToDouble/@StringCases[directory,NumberString]

convertTimeEstimate[time_,ratioIterations_: 10000000/1000]:=ratioIterations*UnitConvert[Quantity[time,"Seconds"],"Hours"]


(* Physical Quantities.  Fundamental constants, relevant length scales. *)

(* Define invariable physical constants: *)
Na=Quantity[6.022*10^23,1/("Moles")];
\[Epsilon]0=Quantity[8.85*10^-12,("Farads")/("Meters")];
\[Epsilon]D=Quantity[78.54*8.85*10^-12,("Farads")/("Meters")];
e=Quantity[1.60217662*10^-19,"Coulombs"];
kB=Quantity[1.3806485279*10^-23,("Kilograms"*("Meters")^2)/(("Seconds")^2*"Kelvins")];
R=Quantity[1.9872*10^-3,("KilocaloriesThermochemical")/("Moles"*"Kelvins")];

(* Bjerrum length, Debye length (and inverse), mean interparticle spacing vs. concentration: *)
\[Lambda]B=UnitConvert[e^2/(4 \[Pi] \[Epsilon]D kB Quantity[298,"Kelvins"]),"Nanometers"];
\[Lambda]D[\[Mu]_]=(1/(8\[Pi] \[Lambda]B Na Quantity[\[Mu],("Moles")/("Decimeters")^3]))^(1/2); 
\[Kappa][\[Mu]_]=\[Lambda]D[\[Mu]]^-1;
spacing[\[Mu]_]=UnitConvert[CubeRoot[(Quantity[\[Mu],"Molar"]Na)^-1],"Nanometers"];

(* Assess the unit of time in LJ units with a given \[Sigma], m, & \[Epsilon]: *)
\[Tau]LJ[\[Sigma]_:56,m_:2*10^7,\[Epsilon]LJ_: (kB Quantity[298,"Kelvins"])]:=Quantity[\[Sigma],"Nanometers"]*Sqrt[Quantity[m,"Grams"/"Moles"]/(\[Epsilon]LJ*Na)]

(* Similarly assess the unit of charge (reduced charge) in LJ units (the following are equivalent, latter overwrites former and is faster: *)
qLJ[z_, \[Sigma]_:56, \[Epsilon]LJ_:(kB Quantity[298,"Kelvins"])]:=UnitSimplify[(z e)/Sqrt[4\[Pi] \[Epsilon]0 Quantity[\[Sigma],"Nanometers"] \[Epsilon]LJ]]

(* An (at least approximate) equation for the net screened ES-energy of a sphere of uniform net charge (from Vanderschoot paper): *)
netScrSphPot[qStr_,cS_,r_:10]:=qStr^2/(4 Quantity[r,"Nanometers"]^2) (\[Lambda]B*\[Lambda]D[cS])

(* An equation for the net surface tension energy of a sphere (in kB T): *)
netTensionPot[\[Sigma]_,r_:10]:=(Quantity[\[Sigma],("Dynes")/("Centimeters")]*(4\[Pi] Quantity[r,"Nanometers"]^2))/(kB*Quantity[298,"Kelvins"])

(* Assesses the number of particles in a given volume (units must be specified) at a given concentration (in Molar): *)
numParticles[conc_,vol_]:=Na*Quantity[conc,"Molar"]*vol

(* The moduli associated with a shell (2D Young's and 3D Bulk modulus, depending on a hollow shell thickness t): *)
youngsModulusY[\[Kappa]s_,R_]:=UnitConvert[(\[Kappa]s*kB*Quantity[298,"Kelvins"])/Quantity[R,"Nanometers"]^2,("Newtons")/("Meters")]
bulkModulusE[\[Kappa]s_,R_,t_]:=UnitConvert[youngsModulusY[\[Kappa]s,R]/Quantity[t,"Nanometers"],"Kilopascals"]


(* Preliminary potential component definitions:  *)

(* The DLVO effective charge as a function of particle diameter (\[Sigma], nanometers), bare charge (z, dimensionless), and ionic strength (\[Mu], Molar):: *)
qDLVO[\[Sigma]_,z_,\[Mu]_]=z(E^(\[Kappa][\[Mu]]*Quantity[(\[Sigma]/2),"Nanometers"])/(1+\[Kappa][\[Mu]]*Quantity[(\[Sigma]/2),"Nanometers"]));

	(* An equivalent dimensionless version (verified 2018.02.07): *)
	(*qDLVO[\[Sigma]_,z_,\[Mu]_]=z(E^((3.2879785861623345/Sqrt[1/\[Mu]])*(\[Sigma]/2))/(1+(3.2879785861623345/Sqrt[1/\[Mu]])*(\[Sigma]/2)));*)


(* Physical Potential Functions: *)

(* The 'hard-core diameter' and hard-core LJ potential: *)
	(* The 'hardcore diameter:'  *)
\[CapitalDelta][\[Sigma]1_,\[Sigma]2_,\[Sigma]HC_]:=(\[Sigma]1+\[Sigma]2)/2-\[Sigma]HC;
	(* The potential itself: *)
uLJ[r_,\[Sigma]1_,\[Sigma]2_,\[Sigma]HC_,\[Epsilon]_:1.]:=Piecewise[{
{\[Infinity],r<=\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]},
{\[Epsilon]+4 \[Epsilon] ((\[Sigma]HC/(r-\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]))^12-(\[Sigma]HC/(r-\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]))^6),\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]<=r<2^(1/6) \[Sigma]HC+\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]},
{0,r>=2^(1/6) \[Sigma]HC+\[CapitalDelta][\[Sigma]1,\[Sigma]2,\[Sigma]HC]}}]

(* The standard, Coulombic potential with no screening (u(1,1,1) = |Subscript[\[Lambda], B]|): *)
uCoul[r_,z1_,z2_]=\[Lambda]B((z1 z2)/Quantity[r,"Nanometers"]);

	(* A cut-off version of the Coulombic potential: *)
	uCoulCut[r_,z1_,z2_,rCut_]:=If[r<rCut,uCoul[r,z1,z2],0.];

(* The screened Yukawa potential: *)
uYuk[r_,z1_,z2_,\[Mu]_]=((\[Lambda]B z1 z2)/Quantity[r,"Nanometers"])*E^(-\[Kappa][\[Mu]] Quantity[r,"Nanometers"]);
	
	(* A cut-off version of the Yukawa potential: *)
	uYukCut[r_,z1_,z2_,\[Mu]_,rCut_]:=If[r<rCut,uYuk[r,z1,z2,\[Mu]],0.];

(* The screened Yukawa potential (with surface shifted DLVO charge): *)
uYukD[r_,\[Sigma]1_,z1_,\[Sigma]2_,z2_,\[Mu]_]=((\[Lambda]B*qDLVO[\[Sigma]1,z1,\[Mu]]*qDLVO[\[Sigma]2,z2,\[Mu]])/Quantity[r,"Nanometers"])*E^(-\[Kappa][\[Mu]]*Quantity[r,"Nanometers"]);

	(* A cut-off version of the Yukawa potential (with surface shifted DLVO charge): *)
	uYukDCut[r_,\[Sigma]1_,z1_,\[Sigma]2_,z2_,\[Mu]_,rCut_]:=If[r<rCut,uYukD[r,\[Sigma]1,z1,\[Sigma]2,z2,\[Mu]],0.];


(* Time-series & analysis utilities: *)

(* A function to produce a dataset weighted by the x-axis intervals between points: *)
weightedDiffData[dataSeries_]:=WeightedData[Last/@dataSeries,Flatten@{1,Differences[First/@dataSeries]}]

(* Function to take the ordinate (y) mean of the specified percentage of the end of the series: *)
endYAvg[series_,percent_:10]:=Module[{percentageSteps=Round[Length@series/percent]},Mean[series[[-percentageSteps;;,2]]]]

(* A function to add a legend tooltip upon mouseover to each curve (uses globally or locally defined 'legends' variable):  *)
toolTipData[data_]:=Table[Tooltip[data[[i]],legends[[i]]],{i,1,Length@legends}]

(* Give the Tally of a list, but report the fraction of elements corresponding to each unique form: *)
nTally[list_]:=Map[#/.{lhs_,rhs_}:>{lhs,N@(rhs/Length@list)}&,Tally[list]]

(* Give the first abscissa (e.g. time) at which the ordinate (e.g. y-value of time-series) has reached 95% the maximum value (estimate of when, say, hyperbola plateaus): *)
detectInitPlateau[series_List]:=First@First@Cases[series,{x_,y_}/;y>=0.95*Max[Last/@series]]


(* Geometric and mathematical utilities: *)

(* An operator to determine if quantities are approximately equal (yields True if the smaller is within \[PlusMinus]1% of the larger value): *)
lhs_\[TildeTilde]rhs_:=Module[{vars},vars=DeleteDuplicates[Cases[{lhs,rhs},s_Symbol/;Not[ValueQ[s]],Infinity]];
Chop[First[Quiet[FindMaximum[Abs[lhs-rhs],Evaluate[Sequence@@vars]]]],.01*Max[{lhs,rhs}]]==0]

(* A function to compute the periodic distance in a cubic box of length (boxLength), returns scalar (Mod auto-threads, Norm converts to scalar): *)
periodicDistance[r1_,r2_,boxLength_]:=Norm@Mod[r1-r2,boxLength,-(boxLength/2)];

(* A function to return the periodic distance matrix of a set of coordinates within a given box: *)
periodicDistanceMatrix[pts_,boxLength_: 1.0]:=Map[Norm,Mod[Outer[Subtract,pts,pts,1],boxLength,-boxLength/2.],{-2}]

(* Computes the gyration tensor of a set of coordinates: *)
gyrationTensorAndShape[coords_]:=
Module[{mx,my,mz,x,y,z,xx,yy,zz,xy,xz,yz},
{mx,my,mz}=Mean[coords];
{x,y,z}=Transpose[coords];
xx=(x-mx).(x-mx);
yy=(y-my).(y-my);
zz=(z-mz).(z-mz);
xy=(x-mx).(y-my);
xz=(x-mx).(z-mz);
yz=(y-my).(z-mz);
\[ScriptCapitalT]=({
 {xx, xy, xz},
 {xy, yy, yz},
 {xz, yz, zz}
})/Length[coords];

{\[Lambda]x,\[Lambda]y,\[Lambda]z}=Sqrt/@Chop/@Sort[Eigenvalues[\[ScriptCapitalT]]];
radiusGyration=Sqrt[Total@{\[Lambda]x^2,\[Lambda]y^2,\[Lambda]z^2}]; (* Equivalent to Sqrt[Tr[\[ScriptCapitalT]]] *)
asphericity=\[Lambda]z^2-((\[Lambda]x^2+\[Lambda]y^2)/2);
normAsphericity=asphericity/radiusGyration^2;
acylindricity=(\[Lambda]y^2-\[Lambda]x^2);
normAcylindricity=acylindricity/radiusGyration^2;
anisotropy=(asphericity^2+(3acylindricity^2)/4)/radiusGyration^4;
ratio=(\[Lambda]x/\[Lambda]z);
\[ScriptCapitalT]]

(* Computes the scalar periodic distance from the (Subscript[Overscript[r, \[RightVector]], PoI]) to the nearest entity (say macroion) in latter list argument: *)
distFromCoordList[PoI_, macroCoords_List, boxLengthReal_]:=
Module[{realDistFromMacros=Norm@periodicDistance[PoI,#,boxLengthReal]&/@macroCoords,minDistFromMacro},
minDistFromMacro=First@Sort@realDistFromMacros]
