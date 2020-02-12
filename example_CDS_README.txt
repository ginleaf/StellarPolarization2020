J/ApJ/872/56    R-band linear polarization of Gaia stars    (Panopoulou+, 2019)
================================================================================
Demonstration of magnetic field tomography with starlight polarization toward
a diffuse sightline of the ISM.
    Panopoulou G.V., Tassis K., Skalidis R., Blinov D., Liodakis I.,
    Pavlidou V., Potter S.B., Ramaprakash A.N., Readhead A.C.S., Wehus I.K.
   <Astrophys. J., 872, 56 (2019)>
   =2019ApJ...872...56P
================================================================================
ADC_Keywords: Polarization; Interstellar medium; Stars, distances; Optical
Keywords: ISM: clouds ; ISM: magnetic fields ; techniques: polarimetric

Abstract:
    The availability of large data sets with stellar distance and
    polarization information will enable a tomographic reconstruction of
    the (plane-of-the-sky-projected) interstellar magnetic field in the
    near future. We demonstrate the feasibility of such a decomposition
    within a small region of the diffuse interstellar medium (ISM). We
    combine measurements of starlight (R-band) linear polarization
    obtained using the RoboPol polarimeter with stellar distances from the
    second Gaia data release. The stellar sample is brighter than 17mag in
    the R-band and reaches out to several kiloparsecs from the Sun. HI
    emission spectra reveal the existence of two distinct clouds along the
    line of sight. We decompose the line-of-sight-integrated stellar
    polarizations to obtain the mean polarization properties of the two
    clouds. The two clouds exhibit significant differences in terms of
    column density and polarization properties. Their mean
    plane-of-the-sky magnetic field orientation differs by 60{deg}. We
    show how our tomographic decomposition can be used to constrain our
    estimates of the polarizing efficiency of the clouds as well as the
    frequency dependence of the polarization angle of polarized dust
    emission. We also demonstrate a new method to constrain cloud
    distances based on this decomposition. Our results represent a preview
    of the wealth of information that can be obtained from a tomographic
    map of the ISM magnetic field.

Description:
    We performed polarimetric observations of our sample during 2016,
    2017, and 2018 with the RoboPol polarimeter, which is mounted on the
    1.3m Ritchey-Chretien telescope at the Skinakas Observatory in Crete,
    Greece. The instrument is an imaging polarimeter, which uses two
    half-wave plates and two Wollaston prisms to simultaneously measure
    the relative Stokes parameters q=Q/I and u=U/I (I is the total
    intensity and Q, U are the absolute Stokes parameters). Observations
    were conducted during 13 nights from 2016 May to July, during five
    nights in 2017 July, and during six nights in 2018 August.
    The observing time for science targets was about 66hr in total.

File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl  Records  Explanations
--------------------------------------------------------------------------------
ReadMe          80        .  This file
table2.dat     139      196  Catalog of stellar polarization measurements
--------------------------------------------------------------------------------

See also:
 II/226 : Stellar polarization catalogs agglomeration (Heiles, 2000)
 I/284  : The USNO-B1.0 Catalog (Monet+ 2003)
 I/337  : Gaia DR1 (Gaia Collaboration, 2016)
 I/347  : Distances to 1.33 billion stars in Gaia DR2 (Bailer-Jones+, 2018)
 I/345  : Gaia DR2 (Gaia Collaboration, 2018)
 J/ApJS/2/389    : Observations of O and B stars (Hiltner, 1956)
 J/ApJS/136/463  : Distances and metallicities of HVCs and IVCs (Wakker, 2001)
 J/A+A/384/1050  : Interstellar polarization. VI. (Berdyugin+, 2002)
 J/ApJ/603/584   : Polarimetry toward Musca dark cloud (Pereyra+, 2004)
 J/A+A/462/621   : UBVRI polarisation in NGC 5749 (Vergne+, 2007)
 J/ApJ/728/104   : Optical polarization for 878 Hipparcos stars (Santos+ 2011)
 J/A+A/561/A24   : Polarization at high galactic latitude (Berdyugin+, 2014)
 J/ApJ/807/5     : R-band polarimetry data in Lupus I region (Franco+, 2015)
 J/MNRAS/452/715 : Optical polarization of the Polaris Flare (Panopoulou+, 2015)
 J/A+A/594/A116  : HI4PI spectra and column density maps (HI4PI team+, 2016)
 J/ApJ/838/80    : JHK polarimetry of stars behind bubble N4 (Chen+, 2017)

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units  Label  Explanations
--------------------------------------------------------------------------------
   1- 19 I19    ---    Gaia   Gaia DR2 (I/345) catalog identifier
  21- 32 A12    ---    USNO   USNO-B1 (I/284) catalog identifier
  34- 42 F9.5   deg    RAdeg  [293.9/295.9] Right Ascension (J2000)
  44- 51 F8.5   deg    DEdeg  [71.8/72.5] Declination (J2000)
  53- 60 F8.5   ---    q      [-0.007/0.04] Stokes q value (Q/I)
  62- 68 F7.5   ---  e_q      [0.0009/0.02] The 1{sigma} uncertainty in q
  70- 77 F8.5   ---    u      [-0.04/0.03] Stokes u value (U/I)
  79- 85 F7.5   ---  e_u      [0.0009/0.02] The 1{sigma} uncertainty in u
  87- 93 F7.5   ---    Pl     [0.0002/0.04] Fractional linear polarization
  95-101 F7.5   ---  e_Pl     [0.0009/0.02] The 1{sigma} uncertainty in Pl
 103-109 F7.5   ---    dPl    [0.0001/0.04] Debiased fractional linear
                               polarization
 111-115 F5.1   ---    theta  [-58/72.4] Polarization angle
 117-120 F4.1   ---  e_theta  [1.4/59.5] The 1{sigma} uncertainty in theta
 122-126 I5     pc     rest   [21/10806]? Estimated distance (r_est)
                               from Bailer-Jones+ (2018, I/347) (1)
 128-131 I4     pc   b_rest   [21/9536]? Lower bound of 68% confidence interval
                               in rest
 133-137 I5     pc   B_rest   [21/12338]? Upper bound of 68% confidence interval
                               in rest
     139 I1     ---    Flag   [0/2]? Source flag (2)
--------------------------------------------------------------------------------
Note (1): Blanks indicate "NaN" values.
Note (2): Source flag as follows:
    0 = intrinsically polarized candidate;
    1 = source is in the 1-Cloud region (mainly Low Velocity Cloud (LVC)
         emission region);
    2 = source is in the 2-Cloud region (significant contribution from
         the Intermediate Velocity Cloud (IVC) and LVC region);
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                     Prepared by [AAS], Emmanuelle Perret [CDS] 03-Jun-2019

