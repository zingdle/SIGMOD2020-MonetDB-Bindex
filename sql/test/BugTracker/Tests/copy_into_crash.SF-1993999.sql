
CREATE TABLE PlateX(
plateID bigint NOT NULL,
spRerun int NOT NULL,
mjd int NOT NULL,
plate smallint NOT NULL,
tile smallint NOT NULL,
mapMjd int NOT NULL,
nExp smallint NOT NULL,
tai real NOT NULL,
raBoresight float NOT NULL,
decBoresight float NOT NULL,
taiHMS varchar(64) NOT NULL,
expTime real NOT NULL,
expTimeB1 real NOT NULL,
expTimeB2 real NOT NULL,
expTimeR1 real NOT NULL,
expTimeR2 real NOT NULL,
helioRV real NOT NULL,
ra real NOT NULL,
"dec" real NOT NULL,
cx float NOT NULL,
cy float NOT NULL,
cz float NOT NULL,
htmID bigint NOT NULL,
sn1_0 real NOT NULL,
sn1_1 real NOT NULL,
sn1_2 real NOT NULL,
sn2_0 real NOT NULL,
sn2_1 real NOT NULL,
sn2_2 real NOT NULL,
dateObs varchar(12) NOT NULL,
timeSys varchar(8) NOT NULL,
quality varchar(12) NOT NULL,
name varchar(32) NOT NULL,
program varchar(16) NOT NULL,
version varchar(64) NOT NULL,
observer varchar(64) NOT NULL,
camVer varchar(64) NOT NULL,
spec2DVer varchar(64) NOT NULL,
utilsVer varchar(64) NOT NULL,
spec1DVer varchar(64) NOT NULL,
readVer varchar(64) NOT NULL,
combVer varchar(64) NOT NULL,
extinction_u real NOT NULL,
extinction_g real NOT NULL,
extinction_r real NOT NULL,
extinction_i real NOT NULL,
extinction_z real NOT NULL,
rOffset1 real NOT NULL,
rSigma1 real NOT NULL,
grOff1 real NOT NULL,
grSigma1 real NOT NULL,
rOffset2 real NOT NULL,
rSigma2 real NOT NULL,
grOff2 real NOT NULL,
grSigma2 real NOT NULL,
sfd_used tinyint NOT NULL,
xygrSig1 real NOT NULL,
xygrSig2 real NOT NULL,
mpgrSig1 real NOT NULL,
mpgrSig2 real NOT NULL,
mpgrOff1 real NOT NULL,
mpgrOff2 real NOT NULL,
isPrimary tinyint NOT NULL,
cartridgeID smallint NOT NULL,
plateVersion varchar(32) NOT NULL,
haMin real NOT NULL,
haMax real NOT NULL,
mjdDesign int NOT NULL,
theta real NOT NULL,
fscanVersion varchar(32) NOT NULL,
fmapVersion varchar(32) NOT NULL,
fscanMode varchar(32) NOT NULL,
fscanSpeed int NOT NULL,
programType int NOT NULL,
programName varchar(32) NOT NULL,
loadVersion int NOT NULL,
expID blob NULL,
CONSTRAINT pk_PlateX_plateID PRIMARY KEY
(
plateID
)
);


COPY 1 RECORDS INTO PlateX from STDIN USING DELIMITERS E'\t', E'\n';
238071330002436096	25	52381	845	616	5237901	12	4.5257452E+9	188.11953	4.9431539999999998	07:03:56.30	3124.0	3.0	3.0	3.0	3.0	10.410612	188.1192	4.9424138	-0.98629544768999999	-0.14070785053962101	8.6154458169632603E-2	14906400923706	19.139601	18.403	15.7439	20.577999	20.612301	17.4648	2000-09-01	tai     	excellent	0845-52379-01	binning 1   1 	v3_74_0 	harvanek	SPEC1 v4_7	v5_3_3  	v5_2_0  	v5_10_4 	v5_3_3  	v5_3_2  	0.1138	8.3800003E-2	6.0699999E-2	4.6100002E-2	3.2699998E-2	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	0	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	-9999.0	1	6	v3_6	6.0	6.0	52336	0.0	0	v4_1_0	extreme	500	0	chunk32	28	4664C8005664D8004664CC005664D8004664D0005664D800580000005664D800580000005664D800580000005664D800578000005664D800578000005664D800578000005664D800584000005664D800584000005664D800584000005664D800


drop table PlateX;
