# Plant Inflorescence counts in the ITEX plots at Barrow and Atqasuk, Alaska

## CONTENT
The data presented are seasonal flowering (number of inflorescences in flower within a plot) collected weekly for all plant species during the summers of 1994-2019 for 48 plots (24 experiment open-top chamber plots and 24 control plots) at four sites (Atqasuk Wet Meadow, Atqasuk Dry Heath, Barrow Wet Meadow, and Barrow Dry Heath).
Please note that raw data for ongoing collections for these datasets (archived in a slightly different format) can be found at: 

Robert Hollister and Katlyn Betway-May. 2022. Plant phenology and performance in the International Tundra Experiment (ITEX) plots at Utqiaġvik and Atqasuk, Alaska. Arctic Data Center. urn:uuid:59f94e1f-93b2-4f8f-b68b-cbf260cd61fc.

## GENERAL INFORMATION

PI/DATA CONTACT= 	Webber, Patrick J (MSU 1994-2005)  
 				Hollister, Robert D (GVSU 2005 onward)    
FUNDING SOURCE AWARD # = NSF 9714103, 0632263, 0856516, 1432277, 1504224, 1836839  
PLATFORM/SITE = Barrow, Alaska (71o31'N 156o60'W);  
                Atqasuk, Alaska (70o45'N 157o40'W)  
INSTRUMENT = Campbell CR10X (MRC TP101M Temperature Probe), Onset HOBO H8 Pro, Onset StowAway XTI, and Onset HOBO XT  

## CONTACT INFORMATION

Robert D. Hollister  
Biology Department  
Grand Valley State University  
1 Campus Drive  
Allendale, MI49401   
Voice: 616 331-8582  
Fax: 616 331-3446  
Email: hollistr [at] gvsu.edu  


## DATA COLLECTION
Plant development was followed throughout the entire summer.  Plant measures were determined based on species morphology and ease of information collection.


### Caveats
Different data were collected in different years.  The data were quality checked, but the vast volume ensures that there are clearly some remaining errors.


## DATA DETAILS

A listing of each species is presented. (these can also be found in tabular form in the file taxon.csv)  

strGSpp	strSpeciesName  
AALP		Alopecurus alpinus  
ABOR		Artemisia borealis  
AFRI		Antennaria friesiana  
AFUL		Arctophila fulva  
ALAT		Arctagrostis latifolia  
BNAN		Betula nana  
CALX		Calamagrostis sp.  
CAQU		Carex aquatilis  
CBEE		Cerastium beeringianum  
CBIG		Carex bigelowii  
CHOL		Calamagrostis holmii  
CHRT		Chrysosplenium tetrandrum  
COFF		Cochlearia officinalis  
CPRA		Cardamine pratensis  
CRAR		Carex rariflora  
CROT		Carex rotundata  
CSTA		Carex aquatilis/stans  
CSUB		Carex subspathacea  
CTET		Cassiope tetragona  
DFIS		Dupontia fisheri  
DLAC		Draba lactea  
DLAP		Diapensia lapponica  
DMIC		Draba micropetala  
DPSI		Dupontia fisheri/psilosantha  
EANG		Eriophorum angustifolium  
ERUS		Eriophorum russeolum  
ESCH		Eriophorum scheuchzeri  
ETRI		Eriophorum angustifolium/triste  
FBRA		Festuca brachyphylla  
HALP		Hierochloe alpina  
HPAU		Hierochloe pauciflora  
JBIG		Juncus biglumis  
LARC		Luzula arctica  
LCON		Luzula confusa  
LPAL		Ledum palustre  
LWAH		Luzula wahlenbergii  
MAPE		Melandrium apetalum  
MOBT		Minuartia obtusiloba  
MUSH		Fungi  
ODIG		Oxyria digyna  

## DATA FORMAT
Headers for the file can be found in tblInfloCountsInfl_headers.txt  
The following prefixes have been used to designate data type.  
   str - indicates string (text) data  
   num - indicates numeric data  

## COLUMN DESCRIPTIONS
strSitCom  - the site and community  
site 1(B=Barrow Alaska, A=Atqasuk Alaska)
community 1(D=Dry Heath Tundra, W=Wet Meadow Tundra)  
strTrea  - treatment 1(C=Control, E=Open-Top Chamber)  
numYear  - year 4(1994-2019)  
strIdLo  - the location of the recording  
	site 1(B=Barrow Alaska, A=Atqasuk Alaska)  
	community 1(D=Dry Heath Tundra, W=Wet Meadow Tundra)  
	treatment 1(C=Control, E=Open-Top Chamber)  
	plot 2(01-24)   
strGSpp  - genus(1) species(3) (See Listing Above)  
strSex  - sex 1(B=Bisexual, M=Male, F=Female, X=Unknown)  
strGSppSex  - genus(1) species(3) (See Listing Above)  
sex 1(B=Bisexual, M=Male, F=Female, X=Unknown)  
strTreaYr  - treatment 1(C=Control, E=Open-Top Chamber)  
year 4(1994-2019)  
numJulian  - day of the year (1-366)  
strJulian  - day of the year (001-366)  
numResponse  - number of inflorescences in flower or withered counted minus the number of withered inflorescences at the previous count (this is the number of inflorescences that flowered since that previous flower count in each plot).  The plots are approximately 1 meter square.

NOTE
       -999.9 - indicates missing data
