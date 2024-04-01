! NUCLEATION ALGORITHM
! AUTHORS: JOSEPH GREGORY Z. CABINTA, ROOSEVELT T. TABAG JR. 
! DATE WRITTEN: MARCH 13, 2019
! DATE UPDATED: OCTOBER 28, 2022
!
! ########################################################################################################
!
!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#      NUCLEATION ALGORITHM        #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#      jgcabinta.rtabag.2022       #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################" 
!
! This algorithm is composed of two main parts.
! First, detection of clusters in an algorithm-readable trajectory file.
! Second, analysis of cluster attributes (behavior and properties).
!
! ########################################################################################################
!
! [ TERMINOLOGIES ] 
!
! star                 -- current molecule serving as the point of reference.~ All  
!                         other molecules are compared to star in the search for clustered molecules.
!                     
! satellite            -- molecules in direct contact with star.~ Remember that in a cluster,
!                         two molecules may not "see" each other and still manage to be in the same
!                         cluster because of the molecules between them. After all molecules are
!                         tested against star, point of reference is transferred to a satellite.
!                         A satellite is called a star satellite when serving as the point of reference.
!                         The term "second-generation satellites" may be used to refer to satellites of star satellite.
!                     
! tolerance            -- vicinity of a star. In the search procedure, a search box is built around the star.~
!                         Only molecules inside the search box is evaluated against Stillinger distance and
!                         energy criteria. 
!
! frag                 -- a unique cluster.~ The total frag in a system is not equal to the total 
!                         number of cluster sizes. A cluster size can have multiple occurrences.
!                         (e.g. system has 3 dimers, 10 trimers, 6 30-mers. The system has
!                         3 cluster sizes and 19 frags).
!                     
! cut molecule/cluster -- molecules or clusters that passed through the plane surface of the simulation box during 
!                         production run.~ Molecular dynamics employs a scheme called periodic boundary condition (PBC). 
!                         PBC is used to approximate a large (infinite) system from a small part called a unit cell.
!                         Replicates of the unit cell is positioned at all of its sides. So when atoms of a molecule 
!                         protrudes the north plane, for example, an identical copy of those atoms appear coming out of
!                         the south plane. This is beneficial to avoid boundary effects in the simulation. However, during
!                         the analysis stage, restoration is necessary to prevent the algorithm from overcounting and 
!                         making the wrong calculations.  
! 
! zero frag            -- frags without any molecule associated to it.~ The process of restoration requires reorganizing 
!                         the unique cluster i.d. assigned to each cluster/frag. Before reorganizing takes place, frags
!                         that will have no molecules associated to it will be present. 
!
! alg                  -- algorithm
!
! ########################################################################################################
!
! [ DESCRIPTOR TAGS AND COUNTERS ] 
!
! Logical Arrays (=.TRUE. or =.FALSE.)
!
! pool(_)             -- TRUE for molecules that are not yet part of a cluster.
!                        This tag is necessary since the search procedure 
!                        doesn't iterate in an organized manner, rather by branching.
!                        So after star matches with a satellite, this tag prevents the 
!                        satellite (when it becomes a star) from scanning star (previous star) 
!                        and repeating the entire clustering evaluation.
!                    
! RNDTAG(_)           -- TRUE for molecules that are evaluated to be part of current cluster.
!                        The final cluster size where star molecule belongs to can only be 
!                        determined after several matchings with other molecules. Hence we cannot 
!                        immediately associate the clustered molecules with a cluster size. So we need 
!                        RNDTAG to tag these molecules first, then apply cluster size and register 
!                        cluster size tally once the entire cluster has been detected. 
!                    
! clus(_)             -- TRUE for satellites that have not yet been made a star.
!                        Each time a satellite has been made a star, clus is simply turned off.
!                        It acts as a list that tells us to keep repeating the search procedure
!                        until no more satellites are left.                                      
!                    
! vapor(_)            -- TRUE for stars that are evaluated to be monomers.
!                        If star is not clustered, it is not involved in the phase change
!                        and remains in the monomer state
!
! Integer Arrays (=N)
!
! nuc(_) = frag       -- This serves as the unique cluster i.d. (sometimes referred to as UCI or frag)
!                        Take note that number of cluster sizes is not equal to total frag
!                        because it is possible to have multiple occurrences of one cluster size.
!                        (e.g. nuc(251)=10; molecule 251 is part of cluster 10)
!
! mer(_) = mercnt     -- cluster size that a molecule is part of. Take note that a monomer is counted as a 1-mer "cluster".  
!                        (e.g. mer(251)=10; molecule 251 is part of a 10-mer) 
!
! tally(_)            -- counts the number of molecules that have a cluster size of $mercnt-mer
!                        (e.g. tally(91)=182; there are 182 molecules that are 91-mers; you can then 
!                        deduce that there are 2 91-mers)  
!
! ########################################################################################################
!
! [ pool, RNDTAG, and clus ] 
! 
! These three descriptor tags may sound very similar. 
! pool is FALSE when molecule is clustered.
! RNDTAG is TRUE when molecule is clustered.
! clus is TRUE when molecule is clustered.
! 
! But their functions are different.
! pool prevents redundancy.
! RNDTAG prepares for uniform attribute tagging.
! clus decides branch reiteration.
!
! pool tag of a molecule will remain FALSE until all molecules in the system has been analyzed.
! RNDTAG tag of a molecule will remain TRUE until all molecules in a cluster has been detected.
! clus tag of a molecule will remain TRUE until molecule (satellite) has finished being a star.
!
! duration of tag: pool > RNDTAG > clus (pool having the longest duration). 
!
! ########################################################################################################
!
! [ HOW TO USE ] 
!
! 1. Set desired parameter settings in fort.10
!
! 2. Clean trajectory file
!          ./ANALYSIS_CLEAN.sh (trajectory filename .w/o extension)
!
! 3. Compile algorithm
!          gfortran -o3 cluster_nm.f -o cluster
!
! 4. Begin analysis
!          ./cluster < (cleaned trajectory file) 
! 
! ########################################################################################################
!
! [ SECTIONS ] 
!
! I. SETTING UP AND REGISTERING NECESSARY PARAMETERS
! I. A. READ PARAMETER FILE
! I. B. INDEX MOLECULES
! I. C. RESTORE CUT CLUSTERS
! I. D. REGISTER CENTER OF GRAVITY
! I. E. ASSIGN DESCRIPTOR TAGS  
! 
! II. BEGIN STILLINGER TEST
! II. A. LOOP ALL SYSTEM MOLECULES
! II. B. GENERATE SEARCH BOX AROUND STAR MOLECULE
! II. C. SCAN SEARCH BOX 
! II. D. LOOP ALL MOLECULES IN POOL (1 to count1)
! II. E. TEST NONBONDED PAIRWISE INTERACTIONS
!
! III. BEGIN BRANCHING 
! III. A. COUNT SATELLITES
! III. B. GENERATE SEARCH BOX AROUND STAR SATELLITE
! III. C. SCAN SEARCH BOX
! III. D. LOOP ALL MOLECULES IN POOL (1 to countb1)
! III. E. TEST NONBONDED PAIRWISE INTERACTIONS
! III. F. DECIDING WHETHER TO FIND NEW STAR OR REITERATE BRANCHING
! III. G. TOP BLOCK TERMINATOR 
!
! IV. BEGIN CLUSTER RESTORATION
! IV. A. COUNT RESTORED MOLECULES
! IV. B. TRANSLATE RESTORED MOLECULES
! IV. C. REGISTER CENTER OF GRAVITY
!
! V. BEGIN STILLINGER TEST FOR RESTORATION 
! V. A. GENERATE SEARCH BOX AROUND STAR MOLECULE 
! V. B. SCAN SEARCH BOX
! V. C. LOOP ALL MOLECULES IN POOL (1 to countp1)
! V. D. TEST NONBONDED PAIRWISE INTERACTIONS
! V. E. RECOVERING CUT CLUSTER
! V. F. REORGANIZING FRAGS
!
! VI. BEGIN PROPERTY AND BEHAVIOR ANALYSIS
!
! VI. A. 1. SEARCHING FOR CUT CLUSTERS
! VI. A. 2. RECOVERING CUT CLUSTER
! 
! VI. B. RADIAL DISTRIBUTION PROFILE #1
! VI. B. 1. MATCHING TARGET TIME
! VI. B. 2. PREPARING REPORT FILES
! VI. B. 3. SETTING STORAGE VARIABLES
! VI. B. 4. SEARCHING CENTER r=0
! VI. B. 5. SEARCHING MOLEC DISTANCE FROM CENTER 
! VI. B. 7. STORING DIST BET CELL CENTER AND INDIV MOLEC
! VI. B. 8. REMOVING EXTRA VALUES
! VI. B. 9. CALCULATING VOLUME AND DENSITY
!
! VI. C.  RADIAL DISTRIBUTION PROFILE #2  
! VI. C.  1. MATCHING TARGET TIME AND MER
! VI. C.  2. PREPARING REPORT FILES
! VI. C.  3. SETTING STORAGE VARIABLES
! VI. C.  4. ITERATE ALL MER SIZE
! VI. C.  5. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. C.  6. INDEXING MOLECULES INVOLVED IN LOCATED FRAG
! VI. C.  7. SEARCHING CENTER OF WHOLE CLUSTER
! VI. C.  8. SEARCHING REFERENCE ATOM OF INDIV MOLEC
! VI. C.  9. DIST BET CLUSTER CENTER AND INDIV MOLEC
! VI. C. 10. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. C. 11. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. C. 12. SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
! VI. C. 13. REMOVING EXTRA VALUES
! VI. C. 14. AVERAGING COUNT TO REPRESENT CHOSEN RANGE
! VI. C. 15. CALCULATING VOLUME AND DENSITY
!
! VI. D.  RADIAL DISTRIBUTION PROFILE #3
! VI. D.  1. MATCHING TARGET TIME AND MER
! VI. D.  2. PREPARING REPORT FILES
! VI. D.  3. SETTING STORAGE VARIABLES
! VI. D.  4. ITERATE ALL MER SIZE
! VI. D.  5. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. D.  6. INDEXING MOLECULES INVOLVED IN LOCATED FRAGS
! VI. D.  7. SEARCHING CENTER OF WHOLE CLUSTER
! VI. D.  8. SEARCHING REFERENCE ATOM OF INDIV MOLEC
! VI. D.  9. DIST BET CLUSTER CENTER AND INDIV MOLEC
! VI. D. 10. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. D. 11. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. D. 12. SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
! VI. D. 13. REMOVING EXTRA VALUES
! VI. D. 14. AVERAGING COUNT TO REPRESENT CHOSEN RANGE
! VI. D. 15. CALCULATING VOLUME AND DENSITY
!
! VI. E.  RADIAL DISTRIBUTION FUNCTION
! VI. E.  1. MATCHING TARGET TIME AND MER
! VI. E.  2. PREPARING REPORT FILES
! VI. E.  3. SETTING STORAGE VARIABLES
! VI. E.  4. ITERATE ALL MER SIZE
! VI. E.  5. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. E.  6. INDEXING MOLECULES INVOLVED IN LOCATED FRAGS
! VI. E.  7. SEARCHING CENTER OF WHOLE CLUSTER
! VI. E.  8. SEARCHING CENTER OF INDIV MOLEC
! VI. E.  9. DIST BET CLUSTER CENTER AND INDIV MOLEC
! VI. E. 10. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. E. 11. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. E. 12. SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
! VI. E. 13. REMOVING EXTRA VALUES
! VI. E. 14. AVERAGING COUNT TO REPRESENT CHOSEN RANGE
! VI. E. 15. CALCULATING VOLUME AND DENSITY
!
! VI. F. CLUSTER SIZE FREQUENCY
! VI. F. 1. SETTING UP VARIABLES
! VI. F. 2. ITERATE ALL MER SIZE
! VI. F. 3. APPLYING PLACE HOLDERS
! VI. F. 4. PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP
!
! VI. G.  MOLE FRACTION  
! VI. G.  1. SETTING UP VARIABLES 
! VI. G.  2. ITERATE ALL MER SIZE
! VI. G.  3. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. G.  4. INDEXING MOLECULES INVOLVED IN LOCATED FRAG
! VI. G.  5. CALCULATING MOLE FRACTION OF CURRENT FRAG
! VI. G.  6. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. G.  7. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. G.  8. SUMMING COUNT TO REPRESENT CHOSEN RANGE
! VI. G.  9. APPLYING PLACE HOLDERS
! VI. G. 10. PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP  
!
! VI. H.  CLUSTER MASS
! VI. H.  1. SETTING UP VARIABLES
! VI. H.  2. ITERATE ALL MER SIZE
! VI. H.  3. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. H.  4. INDEXING MOLECULES INVOLVED IN LOCATED FRAG
! VI. H.  5. CALCULATING MASS OF CURRENT FRAG
! VI. H.  6. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. H.  7. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. H.  8. SUMMING COUNT TO REPRESENT CHOSEN RANGE
! VI. H.  9. APPLYING PLACE HOLDERS
! VI. H. 10. PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP
!
! VI. I.  CLUSTER RADIUS
! VI. I.  1. SETTING UP VARIABLES
! VI. I.  2. ITERATE ALL MER SIZE
! VI. I.  3. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. I.  4. INDEXING MOLECULES INVOLVED IN LOCATED FRAG
! VI. I.  5. CALCULATING CENTER OF MASS OF CURRENT FRAG
! VI. I.  6. CALCULATING CENTER OF MASS OF INDIV MOLEC
! VI. I.  7. CALCULATING RADIUS FROM MOMENT OF INERTIA
! VI. I.  8. SUMMING COUNT FOR CLUSTERS OF SAME SIZE
! VI. I.  9. AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
! VI. I. 10. SUMMING COUNT TO REPRESENT CHOSEN RANGE
! VI. I. 11. APPLYING PLACE HOLDERS
! VI. I. 12. PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP
!
! VI. J. NUCLEATION RATE
! VI. J. 1. SETTING UP VARIABLES
! VI. J. 2. ITERATE ALL MER SIZE
! VI. J. 3. SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE 
! VI. J. 4. ACCOUNTING FRAG IN THRESHOLD BIN 
!
! VI. K.  CRITICAL CLUSTER SIZE
! VI. K.  1. SETTING UP VARIABLES
! VI. K.  2. ITERATE MER SIZE
! VI. K.  3. SEARCHING FRAGS INVOLVED IN CURRENT SIZE
! VI. K.  4. SEARCHING DECAY AND GROWTH OCCURRENCE
! VI. K.  5. DETERMINE ACCEPTABLE COUNT
! VI. K.  6. TALLY NO CHANGE FOR EACH MER SIZE
! VI. K.  7. TALLY CHANGE COUNT FOR EACH MER SIZE
! VI. K.  8. TIME AVERAGE FOR EACH MER SIZE
! VI. K.  9. TRANSITION PROBABILITY AND CHANGE RATE
! VI. K. 10. RECORD MER, FRAG, AND TALLY ID FOR NEXT TIME STEP
!
! VI. L. CLUSTER SNAPSHOTS
! VI. L. 1. MATCHING TARGET TIME AND MER
! VI. L. 2. WRITE TOTAL ATOMS 
! VI. L. 3. ITERATE ALL MER SIZE
! VI. L. 4. SEARCHING FRAGS INVOLVED IN CURRENT SIZE
! VI. L. 5. INDEXING MOLECULES INVOLVED IN LOCATED FRAG
! VI. L. 6. COMMENCE WRITING STRUCTURE FILE
!
! VI. M. MOLECULE REPORT - WRITE MOLECULE INFORMATION
! VI. N. NUCLEI REPORT - WRITE NUCLEI INFORMATION
! VI. O. MER SIZE REPORT - SUMMARY OF MOLECULE FOR EACH CLUSTER SIZE

! ########################################################################################################

	program cluster

! [ USER-DEFINED VARIABLES ] 

	implicit none
	integer i,j,k,l,h,v,bb,dd,eek,mq,cp,cq 
	integer tt,kk,li,lj,nj,nn,cc,bi,bj,hh,fr,pi,pj,mx,cvc
	integer lp,pb,ww,bn,bbn,countpp,qxq,qxqq,my,jj,uu,pp,mp
	integer cnon,cwat,cbut,camm,kksg,kksw,kksn,kksb,kksa,kksm
	integer kksc,kkso
	integer ffi,ndx,bnj,bli,blj,abi,pop,mmi,fff,ggg,cth,smest 
	integer api,pnj,pli,plj,ppp,prr,gap,ppi,ppk,ppj,ppf,ffj
	integer frag,ffmmi,ffmmf,pndx,merp,outp,frags,gggg,vvv        
	integer moltyj,moltybj,moltypj,skips,recov,bncnt,rrg,comp
	integer manc,binc,fill,pas,fcou,qlv
	integer runi,runj,runbi,runbj,runpi,runpj,tfrag,rrh
	integer count1,count2,countseed,countb1,countb2,trcnv
	integer counttransfer,countp1,tmst,allfrag,allfrag1
	integer fragp,fragpp,fragppp,fragb,empfg,newfrag,fragbb
	integer typcnt1,typcnt2,typcnt3,typcnt4,typcnt5,typcnt6
	integer typcnt8
	integer audsum2,audsum,mercnt,endtme,countnonmono
	integer clustersize,mltot,gettfrag,monomer,nonmonomer
	integer vall3,rra,rrf,rwr,rne,rbl,raa
	integer qxv,cgen,lrb,psk,pse
	integer cth05,cth10,cth20,cth30,cth40,cth50
	integer cth15,cth25,cth35,cth45,cth60,cth70,cth80,cth90
	integer tc,nth,np,mim,mfm,nin,nfn,mxv,rmc,rmd,snp,nap
	integer smmi1,smmi2,immi1,immi2,lmmi1,lmmi2,emmix
	integer plcholda,bbx,t1,t2,clock_rate,clock_max,prevfrag
	integer RPRD1,RPRD2,RDF,RPSF,RPMOF,RPM,RPR,RPNR,RPS,RPTR,RPCS
	integer RPRD3
	integer growthsum(1:60000),decaysum(1:60000),cclm2(1:5000)
	integer cclm3(1:5000)
	integer le(1:20),crm(1:60000),cee(1:60000),crm3(1:60000)
	integer crm4(1:60000)
	integer rwat(1:60000),rnon(1:60000),rbut(1:60000),ramm(1:60000)
	integer rmet(1:60000),ract(1:60000),roct(1:60000)
	integer molty(1:60000),ndxx(1:60000),sndx(1:60000),mer(1:60000)
	integer tally(1:60000),nuc(1:60000),pct(1:60000),pndxx(1:60000)
	integer mers(1:60000),val2(1:60000),cclm(1:40)
	integer sumwat(1:20,1:60000),sumnon(1:20,1:60000)
	integer sumbut(1:20,1:60000),sumamm(1:20,1:60000)
	integer summet(1:20,1:60000),sumact(1:20,1:60000)
	integer sumoct(1:20,1:60000)
	integer genmol(1:20,1:60000),kin(1:5000)
	integer kin2(1:5000)
	integer sumwat2(1:20,1:60000),sumnon2(1:20,1:60000)
	integer sumbut2(1:20,1:60000),sumamm2(1:20,1:60000)
	integer summet2(1:20,1:60000),sumact2(1:20,1:60000)
	integer sumoct2(1:20,1:60000)
	integer genmol2(1:20,1:60000)
	integer sumwat3(1:20,1:60000),sumnon3(1:20,1:60000)
	integer sumbut3(1:20,1:60000),sumamm3(1:20,1:60000)
	integer summet3(1:20,1:60000),sumact3(1:20,1:60000)
	integer sumoct3(1:20,1:60000)
	integer genmol3(1:20,1:60000)
	integer sumwat4(1:20,1:60000),sumnon4(1:20,1:60000)
	integer sumbut4(1:20,1:60000),sumamm4(1:20,1:60000)
	integer summet4(1:20,1:60000),sumact4(1:20,1:60000)
	integer sumoct4(1:20,1:60000)
	integer genmol4(1:20,1:60000)
	integer rgen(1:500)
	integer prevnuc(1:60000),prevmer(1:60000)
	integer prevtally(1:60000)
	integer growthtally(1:10000,1:200)
	integer decaytally(1:10000,1:200)
	integer zerotally(1:10000)
	integer decaynuc(1:60000),growthnuc(1:60000)
	integer rwat2(1:60000),rnon2(1:60000),rbut2(1:60000)
	integer ramm2(1:60000),rmet2(1:60000),roct2(1:60000)
        integer rgen2(1:500)
	integer ract2(1:60000)
	integer rwat3(1:60000),rnon3(1:60000),rbut3(1:60000)
	integer ramm3(1:60000),rmet3(1:60000),roct3(1:60000)
        integer rgen3(1:500)
	integer ract3(1:60000)
	integer rwat4(1:60000),rnon4(1:60000),rbut4(1:60000)
	integer ramm4(1:60000),rmet4(1:60000),roct4(1:60000)
        integer rgen4(1:500)
	integer ract4(1:60000)
	double precision e,f,g,a
	double precision ca,cb,cd,ce,cf,cg,ee,yy,ts,tss,bq,br
	double precision sux,suy,suz,qqx,rad,rad2,moi,ms,bnr
	double precision xij,yij,zij,mnx,mny,mnz,mxx,mxy,mxz,bnn
	double precision msm,msx,msy,msz,cmx,cmy,cmz,rrx,rry,rrz,der
	double precision tol,tolb,tolp,ene,coul,coul1,xang,yang,zang
	double precision tmoi,trad,mssm,mssx,mssy,mssz,cmmx,cmmy,cmmz
	double precision massl,massh,relmas,SECONDS,larg,lres
	double precision tidv
	double precision sigx,epsx,mlfrcn1,sigij,epsij,sig6,discut,enecut
	double precision disij,cutoff,sybx,sybx1,sybx2,ecwn,ecnb,ecba
	double precision disij2,temp,difx,dify,difz,avrgrd,lenn,bnnn
	double precision bnstp,avrgms,tmas,avrgup,avrglo,getavml1,tstep
	double precision getavml2,getavrad,getavlom,getavupm,ctme
	double precision uprlim,lorlim,upl,lol,tupl,tlol,getavmas
	double precision mlfrcn2,tmlfrcn1,tmlfrcn2,avrgml1,avrgml2
	double precision mlfrcn3,tmlfrcn3,tmlfrcn4,avrgml3,avrgml4
	double precision mlfrcn4,getavml3,getavml4,brad2
	double precision mlfrcn5,getavml5,tmlfrcn5,avrgml5
	double precision mlfrcn6,getavml6,tmlfrcn6,avrgml6
	double precision mlfrcn8,getavml8,tmlfrcn8,avrgml8
	double precision rrb,vall1,vall2,vall4,vall5,vall6,vall7
	double precision vall2W,vall2N,vall2B,vall2A,vall2M,ovdens
	double precision vall2C,vall2O
	double precision coox,cooy,cooz,cddx,cddy,cddz,cox,coy,coz
	double precision rsx,rsy,rsz,cdx,cdy,cdz,divvol,timme
	double precision plcholdb,plcholdc,plcholdd
	double precision avog,qik,qok,coll 
	double precision gt(1:20)
	double precision ndwat(1:20,1:60000),ndnon(1:20,1:60000)
	double precision ndbut(1:20,1:60000),ndamm(1:20,1:60000)
	double precision ndmet(1:20,1:60000),ndact(1:20,1:60000)
	double precision ndoct(1:20,1:60000)
	double precision ndgen(1:20,1:60000)
	double precision sig(1:20,1:20),eps(1:20,1:20),qq(1:20,1:20)
	double precision mss(1:20,1:20)
	double precision cgx(1:60000),cgy(1:60000),cgz(1:60000)
	double precision x(1:60000,1:20,1:20),y(1:60000,1:20,1:20)
	double precision z(1:60000,1:20,1:20),time(1:60000),plne(1:60000)
	double precision val0(1:60000),val1(1:60000),val3(1:60000)
	double precision val4(1:60000),val5(1:60000),val6(1:60000)
	double precision dvwat(1:20,1:60000),dvnon(1:20,1:60000)
	double precision dvbut(1:20,1:60000),dvamm(1:20,1:60000)
	double precision dvmet(1:20,1:60000),dvact(1:20,1:60000)
	double precision dvoct(1:20,1:60000)
	double precision gmin(1:60000,1:20),gmax(1:60000,1:20)
	double precision snap(1:60000,1:20),dvgen(1:20,1:60000)
	double precision avgzty(1:60000),bankzero
	double precision avgdty(1:60000,1:200),bankgrowth
	double precision avggty(1:60000,1:200),bankdecay,changesum
	double precision tgprob(1:60000,1:200),tdprob(1:60000,1:200)
	double precision derate(1:60000),grrate(1:60000)
	double precision sumchangerate(1:60000)
	double precision rmmw(1:60000),rmmn(1:60000),rmmb(1:60000)
	double precision rmmg(1:60000),rmma(1:60000),rmmm(1:60000)
	double precision rmmc(1:60000)
	double precision rmmw2(1:60000),rmmn2(1:60000),rmmb2(1:60000)
	double precision rmmg2(1:60000),rmma2(1:60000),rmmm2(1:60000)
	double precision rmmc2(1:60000),rmmo2(1:60000)
	double precision rmmw3(1:60000),rmmn3(1:60000),rmmb3(1:60000)
	double precision rmmg3(1:60000),rmma3(1:60000),rmmm3(1:60000)
	double precision rmmc3(1:60000),rmmo3(1:60000)
	double precision rmmw4(1:60000),rmmn4(1:60000),rmmb4(1:60000)
	double precision rmmg4(1:60000),rmma4(1:60000),rmmm4(1:60000)
	double precision rmmc4(1:60000),rmmo4(1:60000)
	character(LEN=2) q,d
	character(LEN=3) b,p,o
	character(LEN=10) sh
	character(LEN=14) spce,spcf
	character(LEN=15) dp1h1,dp2h1,dfh5
	character(LEN=16) dfh1
	character(LEN=17) sphe
	character(LEN=18) spha,sphb,sphc,sphd,dp2h2,dp1h5
	character(LEN=18) sphaa,sphbb,sphcc,sphdd
	character(LEN=19) dp1h4,dfh2,dp2h3
	character(LEN=20) spca,dfh3,csh1
	character(LEN=21) sf,sg
	character(LEN=23) spcb,spcc,spcd
	character(LEN=23) spcaa,csc1
	character(LEN=24) dfc1
	character(LEN=25) spc1,dp2c1
	character(LEN=26) spcbb,spccc,spcdd,dp2h4
	character(LEN=27) sph4 
	character(LEN=28) dp1h3,dfh4 
	character(LEN=30) dp1c2,oni
	character(LEN=33) sphf
	character(LEN=34) dp1h2
	character(LEN=38) dp1c1
	logical pool(1:60000)
	logical vapor(1:60000),clus(1:60000)
	logical nullbox,iterate,RNDTAG(1:60000),pbccheck(1:60000)
    	logical pxmore(1:60000),pymore(1:60000),pzmore(1:60000)
    	logical pxless(1:60000),pyless(1:60000),pzless(1:60000)
   	logical nonzero(1:60000),nogap,WATER,NONANE,BUTANOL,AMMONIA
    	logical FRGPOOL(1:60000),REVD(1:60000),SETLIM,LABEL,SNAPD

    	call system_clock (t1, clock_rate, clock_max)

	write(*,*) " "
	write(*,*) "[2.A] (BRANCH TASK: cluster_nm.f)"
	write(*,*) " "

!-------I. SETTING UP AND REGISTERING NECESSARY PARAMETERS

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#         ANALYSIS SET-UP          #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"

!-------I.A. READ PARAMETER FILE

        !  [ PURPOSE ] In this section, 
        !
        !  alg will read the parameters we set in fort.10.
        !  The first block are production-specific settings.
        !  The second block are algorithm-specific settings.
        !  The third block are non-bonded parameters from force field.
        !  The fourth block tells alg which analysis we want to be applied
        !  The fifth block are settings for radial density analysis and snapshots

	read( 10,* )                              
	read( 10,* ) ffmmi,ffmmf                       ! number of frames. 
	read( 10,* ) outp                              ! coordinate output.   
	read( 10,* ) tstep                             ! timestep (ps).
	read( 10,* ) trcnv                             ! trjconv skip (ps).
	read( 10,* ) pop                               ! population.
	read( 10,* ) temp                              ! temperature.
 	read( 10,* ) sybx                              ! system box length.
 	                                             
 	read( 10,* )                                   
 	read( 10,* )                                   
 	read( 10,* ) tol                               ! box size tolerance of star.     
 	read( 10,* ) tolb                              ! box size tolerance for satellite.
 	read( 10,* ) tolp                              ! box size for translation.
 	read( 10,* ) cutoff                            ! Stillinger dist. cutoff.
 	read( 10,* ) ecwn                              ! Stillinger energy cutoff (water-nonane), (ammonia-nonane).
 	read( 10,* ) ecnb                              ! (nonane-butanol).
 	read( 10,* ) ecba                              ! (butanol-ammonia), (butanol-water), (water-ammonia).
 	
 	read( 10,* )
 	read( 10,* )   
 	do nn = 1,43                                   ! trappe non-bonded parameters.  
 	 read( 10,* ) tt ,kk , qqx , sigx , epsx , ms  ! molec,atom,charge,sigma,epsilon,mass.             
 	  qq( tt,kk ) = qqx                              
 	 sig( tt,kk ) = sigx                           
 	 eps( tt,kk ) = epsx    
 	 mss( tt,kk ) = ms        
 	enddo
 	
 	read( 10,* )                              
 	read( 10,* )                                         
 	read( 10,* ) RPRD1                             ! radial density profile #1 (rprd1).     
 	read( 10,* ) RPRD2                             ! radial density profile #2 (rprd2).
 	read( 10,* ) RPRD3                             ! radial density profile #3 (rprd3).
 	read( 10,* ) RDF                               ! radial distribution function (rdf).
 	read( 10,* ) RPSF                              ! size frequency (rpsf).
 	read( 10,* ) RPMOF                             ! mole fraction (rpmof).
 	read( 10,* ) RPM                               ! mass, mass fraction max, mass fraction min.
 	read( 10,* ) RPR                               ! radius (rpr).
 	read( 10,* ) RPNR                              ! nucleation rate (rpnr).
 	read( 10,* ) RPCS                              ! critical cluster size (rpcs).
 	read( 10,* ) RPS                               ! snapshots (rps).
 	read( 10,* ) RPTR                              ! timestep report (rptr).
 	                                             
 	read( 10,* )                                 
 	read( 10,* )                                 
 	do np = 1,5                                    ! number of timesteps in which radial analysis will be initiated.
 	 read( 10,* ) ts , tc                        
 	 gt( np ) = ts                                 ! specific timestep.
 	 le( np ) = tc                                 ! specific cluster size range in target timestep (e.g. three ranges in one timestep).
 	 do nth = 1,le( np )                          
 	  read( 10,* ) mim,mfm,snp                   
 	  gmin( np,nth ) = mim                         ! minimum cluster size.
 	  gmax( np,nth ) = mfm                         ! maximum cluster size.
 	  snap( np,nth ) = snp                         ! specific cluster size within range which will be extracted from trajectory.
 	 enddo
 	 read( 10,* )
 	enddo

!-------I.B. INDEX MOLECULES

        !  [ PURPOSE ] In this section, 
        !
        !  alg begins reading the trajectory file we fed it.
        !  It will go through each molecule $i until it reaches $pop.
        !  For each molecule $i, it starts with the first atom.
        !  From here, it determines the identity of the molecule.
        !  By knowing the identity, alg will know how many atoms
        !  comprise the molecule. Hence, it will know how many more 
        !  rows (identified as $l) should it read before considering one 
        !  molecule finished.
        !
        !  After determining the identity and length, it assigns the coordinates to
        !  an array containing (molecule number, molecule identity, atom number). 
        !  This will make it possible to call on the coordinates of a molecule
        !  in a later time.

	do fr = ffmmi,ffmmf                                    
	 ctme = DBLE(fr)*DBLE(outp)*tstep*(DBLE(trcnv)/1000d0)   ! Prints current timestep being analyzed.
   	 do cc = 1,pop                                       
	  tally( cc ) = 0                                    ! Recall, tally counts the number of molecules that have a cluster size of $mercnt-mer.            
	 enddo                                               
                                                             
	 do i = 1,pop                                        ! Iterates through all molecules.          
	  read( 5,* ) b , e , f , g                          ! Reading trajectory file.                   
	  if ( b .eq. "SOL" ) then                           ! Determine molecule identity.           
	   molty(i) = 1                                                
	   l = 4                                             ! Assign length of molecule.   
	  elseif ( b .eq. "NON" ) then                        
	   molty(i) = 2                         
	   l = 9
	  elseif ( b .eq. "BUT" ) then
	   molty(i) = 3                    
	   l = 6
	  elseif ( b .eq. "NH3" ) then
	   molty(i) = 4                        
	   l = 5
	  elseif ( b .eq. "MET" ) then
	   molty(i) = 5                        
	   l = 3
	  elseif ( b .eq. "ACT" ) then
	   molty(i) = 6                        
	   l = 5
    	  elseif ( b .eq. "ARG" ) then      
    	   molty(i) = 7
    	   l = 1
	  elseif ( b .eq. "OCT" ) then
	   molty(i) = 8                        
	   l = 10
    	  endif
    	  x( i,molty(i),1 ) = e                              ! Assigning x coordinate to a molecule number, molecule identity, atom number (in this case, first atom).
    	  y( i,molty(i),1 ) = f                              ! Assigning y coordinate...
    	  z( i,molty(i),1 ) = g                              ! Assigning z coordinate...
                                                           
    	  if ( l .gt. 1 ) then                               ! For non-argon molecules. 
   	   do k = 2,l                                        ! Read coordinates of second atom to last atom of molecule.                 
	    read( 5,* ) o , xang , yang , zang                      
	    x( i,molty(i),k ) = xang         
   	    y( i,molty(i),k ) = yang
	    z( i,molty(i),k ) = zang
	   enddo
    	  endif
                           
!-------I.C. RESTORE CUT MOLECULES

        !  [ PURPOSE ] 
        !
        !  Molecules can be cut by pbc. To identify molecules that were
        !  cut, we measure the distance of the first atom with all other
        !  atoms of the same molecule. If the distance is great (almost
        !  spanning the box length), then we know it has been cut. The
        !  cut atoms are translated an entire box length to fix this. 

    	  sybx1 = sybx - 10                                  ! I assumed that the distance of the cut won't exceed 10 nm.
    	  if ( l .gt. 1 ) then                               ! For non-argon molecules.
    	   do hh = 2,l
    	    difx = x( i,molty(i),1 ) - x( i,molty(i),hh )    ! Distance between first atom and atom 2... 3... 4... etc in x-axis.
    	    dify = y( i,molty(i),1 ) - y( i,molty(i),hh )
    	    difz = z( i,molty(i),1 ) - z( i,molty(i),hh )
    	    if ( abs( difx ) .gt. sybx1 ) then               ! (e.g.) ATOM_1 =1.2, ATOM_2 =39.4; 1.2-39.4 =-38.2; abs(-38.2) =38.2; hence greater than sybx1 =30 (atom too far away from each other)
    	     if ( x( i,molty(i),hh ) .gt. sybx1 ) then       ! Is 39.4 gt than 30?
    	      x( i,molty(i),hh ) = x( i,molty(i),hh ) -sybx  ! Yes, so we subtract 39.4 by sybx =40; 39.4-40 =-0.6 
   	     else                                            ! Hence atom is translated to -0.6 nm of x axis.
    	      x( i,molty(i),hh ) = x( i,molty(i),hh ) +sybx  
    	     endif
    	    endif
    	    if ( abs( dify ) .gt. sybx1 ) then               ! (e.g.) ATOM_1 =39.4, ATOM_2 =1.2; 39.4-1.2 =38.2; abs(38.2) =38.2; hence greater than sybx1 =30 (atoms too far away from each other)
    	     if ( y( i,molty(i),hh ) .gt. sybx1 ) then       ! Is 1.2 gt than 30?
    	      y( i,molty(i),hh ) = y( i,molty(i),hh ) -sybx  
    	     else                                            ! No, so we add 1.2 by sybx =40; 1.2+40 =41.2
   	      y( i,molty(i),hh ) = y( i,molty(i),hh ) +sybx  ! Hence atom is translated to 41.2 nm of y axis.
    	     endif
    	    endif
    	    if ( abs( difz ) .gt. sybx1 ) then 
    	     if ( z( i,molty(i),hh ) .gt. sybx1 ) then
    	      z( i,molty(i),hh ) = z( i,molty(i),hh ) -sybx
   	     else
    	      z( i,molty(i),hh ) = z( i,molty(i),hh ) +sybx
    	     endif
   	    endif
   	   enddo
 	  endif

   	  pbccheck(i) =.FALSE.              ! Indexing molecules exceeding box.
    	  pxmore(i) = .FALSE.               ! The fix earlier was for cut molecules,
    	  pymore(i) = .FALSE.               ! This tags will serve to restore cut clusters later.
    	  pzmore(i) = .FALSE.               ! I'm using the fact that if one molecule was cut,
    	  pxless(i) = .FALSE.               ! it's likely that an entire cluster was cut as well.
    	  pyless(i) = .FALSE.
    	  pzless(i) = .FALSE.

   	  sybx2 = sybx - 10
   	  if ( l .gt. 1 ) then              ! For non-argon molecules. 
    	   do h = 1,l
    	    if ( x( i,molty(i),h ) .gt. sybx ) then                    
    	     pbccheck(i) = .TRUE.
    	     pxmore(i)   = .TRUE.
    	    endif
    	    if ( y( i,molty(i),h ) .gt. sybx ) then 
    	     pbccheck(i) = .TRUE.
    	     pymore(i)   = .TRUE.
    	    endif
   	    if ( z( i,molty(i),h ) .gt. sybx) then
    	     pbccheck(i) =.TRUE.
    	     pzmore(i)   =.TRUE.
    	    endif
    	    if ( x( i,molty(i),h ) .lt. 0 ) then
    	     pbccheck(i) =.TRUE.
    	     pxless(i)   =.TRUE.
   	    endif
    	    if ( y( i,molty(i),h ) .lt. 0 ) then 
    	     pbccheck(i) =.TRUE.
    	     pyless(i)   =.TRUE.
    	    endif
    	    if ( z( i,molty(i),h ) .lt. 0 ) then
   	     pbccheck(i) =.TRUE.
    	     pzless(i)   =.TRUE.
    	    endif
   	   enddo
   	  else
    	   if ( x( i,molty(i),h ) .gt. sybx2 ) then
    	    pbccheck(i) =.TRUE.
   	    pxmore(i)   =.TRUE.
    	   endif
    	   if ( y( i,molty(i),h ) .gt. sybx2 ) then 
    	    pbccheck(i) =.TRUE.
    	    pymore(i)   =.TRUE.
    	   endif
    	   if ( z( i,molty(i),h ) .gt. sybx2 ) then
    	    pbccheck(i) =.TRUE.
   	    pzmore(i)   =.TRUE.
   	   endif
    	   if ( x( i,molty(i),h ) .lt. 1.5 ) then
   	    pbccheck(i) =.TRUE.
    	    pyless(i)   =.TRUE.
   	   endif
   	   if ( y( i,molty(i),h ) .lt. 1.5 ) then
    	    pbccheck(i) =.TRUE.
   	    pyless(i)   =.TRUE.
    	   endif
  	   if ( z( i,molty(i),h ) .lt. 1.5 ) then
	    pbccheck(i) =.TRUE.
    	    pzless(i)   =.TRUE.
   	   endif
   	  endif

!-------I.D. REGISTER CENTER OF GRAVITY

        !  [ PURPOSE ] 
        ! 
        !  To set the center of gravity, we simply get the average of the
        !  coordinates in each axis. 

	  sux = 0d0                             ! Setting up center of gravity.                          
	  suy = 0d0             
	  suz = 0d0
    	  do h = 1,l
	   sux = sux + x( i,molty(i),h )                 
	   suy = suy + y( i,molty(i),h )
	   suz = suz + z( i,molty(i),h )
    	  enddo
	  cgx(i) = sux/DBLE(l)   
	  cgy(i) = suy/DBLE(l)  
	  cgz(i) = suz/DBLE(l)
                                                ! Tagging molecule descriptors. 
!-------I.E. ASSIGN DESCRIPTOR TAGS
        
        !  [ PURPOSE ]  
        ! 
        !  Several descriptor tags (conditional variables) are employed in the algorithm,
        !  they are used to determine the molecule attributes and allows alg to evaluate
        !  what to do with the molecules depending on their tags  

	  pool(i)   = .TRUE.                             
	  vapor(i)  = .FALSE.                     
	  clus(i)   = .FALSE.                   
	  RNDTAG(i) = .FALSE.                   
	 enddo                                  
                                               
!-------II. BEGIN STILLINGER TEST

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#         STILLINGER TEST          #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"

!-------II.A. LOOP ALL SYSTEM MOLECULES
 
        !  [ PURPOSE ] In this section, 
        !
        !  Alg will begin to go through all the molecules and analyze them with the Stillinger test.
        !  Analysis begins with a star molecule.
        !  Alg then scans vicinity (5 nm) from star. Only those within tolerance 5 nm will be subjected to Stillinger test.
        !  Once all molecules within the vicinity are measured, alg will transfer to a satellite.
        !  Then the entire process repeats (look for molecules within vicinity of satellite... Stillinger test).

 	 frag=0                                                 ! Counter for the number of unique clusters. 
	 do i=1,pop                                             ! Iterates through all molecules.                    
	  mercnt=1                                              ! Lowest cluster size is 1.
	  if ( pool(i) .eqv. .TRUE. ) then   
	   oni='(a15,i6,a15,i5,a15,f10.4,a3)'
	  write(*,oni) " STAR MOLECULE:",i,"FRAME:",fr,"TIME:",ctme,"ns"
	   pool(i)=.FALSE.                                      ! FALSE, molecule is under analysis. It should not be analyzed once alg transfers to a star satellite.
	   RNDTAG(i)=.TRUE.                                     ! TRUE, molecule is part of current cluster.
	   frag=frag+1                                          ! Regardless whether star gets matched with another molecule, it is already a cluster (1-mer). 
	   nuc(i)=frag                                          ! Assign the unique cluster i.d. to star.

!-------II.B. GENERATE SEARCH BOX AROUND STAR MOLECULE

	   mnx=cgx(i)-tol                                       ! Only molecules within tolerance will be evaluated against Stillinger distance & energy criteria. 
	   mny=cgy(i)-tol                   
	   mnz=cgz(i)-tol                       
	   mxx=cgx(i)+tol                       
	   mxy=cgy(i)+tol                       
	   mxz=cgz(i)+tol                       

        !  [ QUESTION ] tol as 1.5 instead of 5?  
        !
        !  Can't we immediately make tol =1.5 (which is the Stillinger distance)?
        !  In measuring Stillinger distance, more calculations are necessary than just subtraction. 
        !  Hence, this approach may indeed be better. 
        !  {Stillinger criteria over all molecules} vs. {subtraction criteria over all molecules and Stillinger over few molecules}

!-------II.C. SCAN SEARCH BOX

	   count1=0                                             ! Counter for number of molecules within tolerance. 
	   do j=1,pop                                           ! Compare star to molecule pool.
	    if ( (pool(j) .eqv. .TRUE.)                         ! If molecule is in pool, it is not yet clustered.
     &      .and. ((cgx(j) .gt. mnx) .and. (cgx(j) .lt. mxx))   ! Check if molecule $j is within tolerance.
     &      .and. ((cgy(j) .gt. mny) .and. (cgy(j) .lt. mxy))
     &      .and. ((cgz(j) .gt. mnz) .and. (cgz(j) .lt. mxz)) ) then 
	     count1=count1+1                                
	     ndxx(count1)=j                                     ! Associating a new i.d. for molecule in search box to original molecule i.d. 
	    endif                                                     
	   enddo                                                   
   	                                                       
   	   if ( count1 .eq. 0 ) then                            ! Alg did not find a match. No molecule found within tolerance.
	    vapor(i)=.TRUE.                                     ! Star is a vapor, not clustered.
	    RNDTAG(i)=.FALSE.                                   ! Reset RNDTAG in preparation for next round of analysis.
	    mer(i)=1                                            ! Associate star with its cluster size, which is 1 (Cluster is a 1-mer).
	    tally(1)=tally(1)+1                                 ! Counting number of molecules that are monomers. 
	    GOTO 400                                            ! No need to evaluate for Stillinger distance & energy. No need to initiate branching.
	   else                              

!-------II.D. LOOP ALL MOLECULES IN POOL (1 to count1)
           
	    count2=0
	    do nj=1,count1                                      ! Iterate through all atoms within tolerance.
	     ndx=ndxx(nj)                                       ! This row and next row is a theme that will reoccur throughout entire alg, it is just associating
	     moltyj=molty(ndx)                                  ! a new i.d. to original molecule i.d. to suit the current analysis

!-------II.E. TEST NONBONDED PAIRWISE INTERACTIONS

	     if ( molty(i) .eq. 1 ) then                        ! Calling identity of molecule $i.
	      li=4                                              ! 1 =SOL; 2 =NON; 3 =BUT; 4 =NH3.
	     elseif ( molty(i) .eq. 2 ) then
	      li=9
	     elseif ( molty(i) .eq. 3 ) then
	      li=6
	     elseif ( molty(i) .eq. 4 ) then
	      li=5
    	     elseif ( molty(i) .eq. 5 ) then
    	      li=3
    	     elseif ( molty(i) .eq. 6 ) then
    	      li=5
    	     elseif ( molty(i) .eq. 7 ) then
    	      li=1
    	     elseif ( molty(i) .eq. 8 ) then
    	      li=10
	     endif

	     if ( moltyj .eq. 1) then                           ! Calling identity of molecule $j.                 
	      lj=4
	     elseif ( moltyj .eq. 2 ) then
	      lj=9
	     elseif ( moltyj .eq. 3 ) then
	      lj=6
	     elseif ( moltyj .eq. 4 ) then
	      lj=5
    	     elseif ( moltyj .eq. 5 ) then
   	      lj=3
    	     elseif ( moltyj .eq. 6 ) then
   	      lj=5
    	     elseif ( moltyj .eq. 7 ) then
   	      lj=1
    	     elseif ( moltyj .eq. 8 ) then
   	      lj=10
	     endif
	     nullbox=.TRUE.              
	     do runi=1,li                                               ! Each atom ($runi) of molecule $i compared with 
	      do runj=1,lj                                              ! each atom ($runj) of molecule $j.
	       sigij=0.5d0*(sig(molty(i),runi)+sig(moltyj,runj))        ! Lorentz-Berthelot combining rules.
	       discut=cutoff*sigij                                      ! Distance cut-off length.
	       xij=x(i,molty(i),runi)-x(ndx,moltyj,runj)              
	       yij=y(i,molty(i),runi)-y(ndx,moltyj,runj)
	       zij=z(i,molty(i),runi)-z(ndx,moltyj,runj)
	       disij2=xij**2+yij**2+zij**2
   	       disij=disij2**0.5d0
	       if ( disij .lt. discut ) then                            ! Checking whether molecule is within Stillinger distance cutoff.
	        nullbox=.FALSE.                                         ! If FALSE, distance test passed.
	        EXIT                                                    ! No need to test other atoms. Proceed to Stillinger energy criterion.              
	       endif
	      enddo                                        

	      if ( nullbox .eqv. .FALSE. ) then           
	       EXIT                                      
	      endif
	     enddo                                                       
	  
	     if ( nullbox .eqv. .TRUE. ) then                           ! Stillinger distance failed. Skip Stillinger energy test (essentially proceeding to next molecule in search box). 
	      GOTO 40                                                   ! Even without this line, a conditional is present before energy test that makes sure only 
	     endif                                                      ! those that passed may enter energy test.
   	
	     if ( disij .lt. discut ) then                              ! Commence Stillinger energy test for molecules that passed distance criterion. 
	       if (((molty(i) .eq. 1) .and. (moltyj .eq. 2))            ! We need to call molecule identity because energy criteria will be different for each 
     &        .or. ((molty(i) .eq. 2) .and. (moltyj .eq. 1))            ! molecule interacting.
     &        .or. ((molty(i) .eq. 2) .and. (moltyj .eq. 4))            ! [ water-nonane ][ nonane-water ]  
     &        .or. ((molty(i) .eq. 4) .and. (moltyj .eq. 2))            ! [ nonane-ammonia ][ ammonia-nonane ]
     &        .or. ((molty(i) .eq. 2) .and. (moltyj .eq. 5))            ! [ nonane-methanol ] 
     &        .or. ((molty(i) .eq. 5) .and. (moltyj .eq. 2))            ! [ methanol-nonane ] 
     &        .or. ((molty(i) .eq. 2) .and. (moltyj .eq. 6))            ! [ nonane-acetic_acid ]
     &        .or. ((molty(i) .eq. 6) .and. (moltyj .eq. 2)))  then     ! [ acetic_acid-nonane ]
	        enecut=ecwn/temp                                        ! $ecwn energy value is used here (recall this was set in fort.10).
	      elseif (((molty(i) .eq. 2) .and. (moltyj .eq. 3))         ! [ nonane-butanol ]
     &           .or. ((molty(i) .eq. 3) .and. (moltyj .eq. 2))         ! [ butanol-nonane ]
     &           .or. ((molty(i) .eq. 2) .and. (moltyj .eq. 8))         ! [ nonane-octanol ]
     &           .or. ((molty(i) .eq. 8) .and. (moltyj .eq. 2))) then   ! [ octanol-nonane ]
    	        enecut=ecnb/temp                                        ! $ecnb energy value is used here.
	      else
	        enecut=ecba/temp                                        ! $ecba energy value is used here.
	      endif

	      ene=0d0
	      do runi=1,li                                              ! Each atom ($runi) of molecule $i compared with
	       do runj=1,lj                                             ! each atom ($runj) of molecule $j.
	        xij=x(i,molty(i),runi)-x(ndx,moltyj,runj) 
	        yij=y(i,molty(i),runi)-y(ndx,moltyj,runj)
	        zij=z(i,molty(i),runi)-z(ndx,moltyj,runj)
	        disij2=xij**2+yij**2+zij**2
	        disij=disij2**0.5d0
	        epsij=(eps(molty(i),runi)*eps(moltyj,runj))**0.5d0      ! Lorentz-Berthelot combining rules.
	        sigij=0.5d0*(sig(molty(i),runi)+sig(moltyj,runj))       ! Lorentz-Berthelot combining rules.
	        sig6=(sigij/disij)**6                                  
	        coul1=qq(molty(i),runi)*qq(moltyj,runj)/disij          
	        coul=138.93545*coul1                                    ! Coulumbic potential.
	        lenn=4*epsij*sig6*(sig6-1d0)                            ! Lennard-Jones potential.
	        ene=ene+lenn+coul                                       ! Total intermolecular interaction potential.
    	       enddo                                                   
    	      enddo                                                    
    	
	      if ( ene .lt. enecut ) then                               ! Alg found a match! Molecule passed criteria.
	       nuc(ndx)=frag                                            ! Molecule $ndx assigned the unique cluster i.d.
	       vapor(ndx)=.FALSE.                                       ! Molecule $ndx not a vapor as well since it has been clustered.
	       pool(ndx)=.FALSE.                                        ! FALSE, because molecule we compared to star has already been analyzed.
	       RNDTAG(ndx)=.TRUE.                                       ! TRUE, so later we know which molecules will star properties be applied to.
	       clus(ndx)=.TRUE.                                         ! Molecule $ndx tagged as satellite. 
	       mercnt=mercnt+1                                          ! We add one to current mer count to register newly found molecule in cluster with star.
	       count2=count2+1                                          ! Total number of satellites. 
	       vapor(i)=.FALSE.                                         ! FALSE, because molecule $i not vapor and is clustered.  
	      endif                  
	     endif                  

   40        CONTINUE                                                   ! End point for Stillinger distance failed. Try another molecule in tolerance.
	    enddo                                                       ! End scan of molecules in tolerance.
   	
	    if ( count2 .eq. 0 ) then                                   ! Alg did not find a match. None passed Stillinger distance & energy criteria
	     vapor(i)=.TRUE.                                            ! Star is a vapor, not clustered.
	     RNDTAG(i)=.FALSE.                                          ! Reset RNDTAG in preparation for next round of analysis.
	     nuc(i)=frag                                                ! Assign the unique identifier (cluster i.d.) to star.                                          
	     mer(i)=1                                                   ! Associate star with its cluster size, which is 1 (Cluster is a 1-mer).
	     tally(1)=tally(1)+1                                        ! Counting number of molecules that are monomers.
	     GOTO 200                                                   ! No need to initiate branching. 
	    else                                                       
	     GOTO 100                                                   ! Initiate branching.
	    endif

        !  [ QUESTION ] RNDTAG(i) vs. clus(i) 
        ! 
        !  It seems that RNDTAG(i) and clus(i) is the same?
        !  RNDTAG will mark molecules that matched with star to know which molecules will star properties be applied to
        !  clus will mark molecules that matched with star to assign as satellite
        !  They are similar, but not the same because clus is called much earlier.
        !  As soon as alg will proceed to use satellite, clus is called.
        !  After no more satellites left, it is only then that RNDTAG is called.
         
  100   CONTINUE                                                 

!-------III. BEGIN BRANCHING

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#    PERFORM  BRANCH  ALGORITHM    #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"

!-------III.A. COUNT SATELLITES

        !  [ PURPOSE ] In this section, 
        !  
        !  Satellites will be compared to the molecule pool
        !  Same procedure, start by determining molecules within vicinity
        !  Then applying Stillinger Test (Distance and Energy) 

	    countseed=0
	    do bi=1,pop                                                 
	     if ( clus(bi) .eqv. .TRUE. ) then                        ! Check if molecule is a satellite.                 
	      countseed=countseed+1                                   ! Counter for number of satellites.
	      sndx(countseed)=bi                                      ! All molecule i.d. of satellites are associated with a satellite i.d.
	     endif
	    enddo
c	    write(*,*) "Count satellites:",countseed

!-------III.B. GENERATE SEARCH BOX AROUND STAR SATELLITE

	    do abi=1,countseed                                        ! Iterate through satellites, $abi is current star satellite. 
	     bi=sndx(abi)                                             ! We use satellite i.d. to call molecule i.d.
   	                                                              ! Important because all the properties are connected to molecule i.d.
	     mnx=cgx(bi)-tolb                                         ! Only molecules within tolerance will be evaluated against Stillinger distance & energy criteria.
   	     mny=cgy(bi)-tolb
  	     mnz=cgz(bi)-tolb
   	     mxx=cgx(bi)+tolb
    	     mxy=cgy(bi)+tolb
    	     mxz=cgz(bi)+tolb

!-------III.C. SCAN SEARCH BOX

   	     countb1=0                                                ! Counter for number of molecules within tolerance. 
   	     do bj=1,pop                                              ! Compare satellite to molecule pool.
              if ( (pool(bj) .eqv. .TRUE.)                            ! If molecule is in pool, it is not yet clustered.
     &   .and. ((cgx(bj) .gt. mnx) .and. (cgx(bj) .lt. mxx))          ! Check if molecule $bj is within tolerance.
     &   .and. ((cgy(bj) .gt. mny) .and. (cgy(bj) .lt. mxy))
     &   .and. ((cgz(bj) .gt. mnz) .and. (cgz(bj) .lt. mxz)) ) then
    	       countb1=countb1+1                             
    	       ndxx(countb1)=bj                                       ! Associating branch molecule i.d. to original molecule i.d.
    	      endif                           
    	     enddo                          
   	
   	     if (countb1 .eq. 0 ) then                                ! Alg did not find a match. No molecule found within tolerance.          
    	      clus(bi)=.FALSE.                                        ! Remove satellite tag on star satellite and move on to a new star satellite.
    	     else                                                     

!-------III.D. LOOP ALL MOLECULES IN POOL (i to countb1)

	      countb2=0                                               ! Counter for number of molecules that will pass Stillinger distance & energy criteria.
	      do bnj=1,countb1                                        ! Iterate through all molecule that passed tolerance. 
	       ndx=ndxx(bnj)                                          ! This row and next row is a theme that will reoccur throughout entire alg, it is just associating 
	       moltybj=molty(ndx)                                     ! a new i.d. to original molecule i.d. to suit the current analysis.

!-------III.E. TEST NONBONDED PAIRWISE INTERACTIONS

	       if ( molty(bi) .eq. 1 ) then                           ! Calling identity of molecule $bi  
	        bli=4                                                 ! 1 =SOL; 2 =NON; 3=BUT; 4 =NH3
	       elseif ( molty(bi) .eq. 2 ) then
	        bli=9
	       elseif ( molty(bi) .eq. 3 ) then
	        bli=6
	       elseif ( molty(bi) .eq. 4 ) then
	        bli=5
    	       elseif ( molty(bi) .eq. 5 ) then
   	        bli=3
    	       elseif ( molty(bi) .eq. 6 ) then
   	        bli=5
    	       elseif ( molty(bi) .eq. 7 ) then
   	        bli=1
    	       elseif ( molty(bi) .eq. 8 ) then
   	        bli=10
	       endif
   	
	       if ( moltybj .eq. 1) then                               ! Calling identity of molecule $bj.
	        blj=4
	       elseif ( moltybj .eq. 2 ) then
	        blj=9
	       elseif ( moltybj .eq. 3 ) then
	        blj=6
	       elseif ( moltybj .eq. 4 ) then
	        blj=5
    	       elseif ( moltybj .eq. 5 ) then
    	        blj=3
    	       elseif ( moltybj .eq. 6 ) then
    	        blj=5
    	       elseif ( moltybj .eq. 7 ) then
    	        blj=1
    	       elseif ( moltybj .eq. 8 ) then
    	        blj=10
	       endif
   	
	       nullbox=.TRUE.
	       do runbi=1,bli                                           ! Each atom ($runbi) of molecule $bi compared with 
	        do runbj=1,blj                                          ! each atom ($runbj) of molecule $bj.
	         sigij=0.5d0*(sig(molty(bi),runbi)+sig(moltybj,runbj))  ! Lorentz-Berthelot combining rules.
	         discut=cutoff*sigij                                    ! Distance cut-off length.
	         xij=x(bi,molty(bi),runbi)-x(ndx,moltybj,runbj)         
	         yij=y(bi,molty(bi),runbi)-y(ndx,moltybj,runbj)         
	         zij=z(bi,molty(bi),runbi)-z(ndx,moltybj,runbj)         
	         disij2=xij**2+yij**2+zij**2                            
	         disij=disij2**0.5d0                                    
	         if ( disij .lt. discut ) then                          ! Checking whether molecule is within Stillinger distance cutoff.
	          nullbox=.FALSE.                                       ! If FALSE, distance test passed. 
	          EXIT                                                  ! No need to test other atoms. Proceed to Stillinger energy criterion. 
	         endif                                                
	        enddo                                                 
   	
	        if ( nullbox .eqv. .FALSE. ) then        
	         EXIT
	        endif
	       enddo                                                 
    	       
	       if ( nullbox .eqv. .TRUE. ) then                         ! Stillinger distance failed, Skip Stillinger energy test (essentially proceeding to next molecule in search box).
	        clus(bi)=.FALSE.                                        ! Turn off satellite tag, finished as star satellite.
	        GOTO 300                                                ! Even without this line, a conditional is present before energy test that makes sure only
	       endif                                                    ! those that passed may enter energy test. 
   	                                                            
	       if ( disij .lt. discut ) then                            ! Commence Stillinger energy test for molecules that passed distance criterion. 
	         if (((molty(bi) .eq. 1) .and. (moltybj .eq. 2))        ! We need to call molecule identity because energy criteria will be different for each 
     &          .or. ((molty(bi) .eq. 2) .and. (moltybj .eq. 1))        ! molecule interacting.
     &          .or. ((molty(bi) .eq. 2) .and. (moltybj .eq. 4))
     &          .or. ((molty(bi) .eq. 4) .and. (moltybj .eq. 2))
     &          .or. ((molty(bi) .eq. 2) .and. (moltybj .eq. 5))        ! [ nonane-methanol ]
     &          .or. ((molty(bi) .eq. 5) .and. (moltybj .eq. 2))        ! [ methanol-nonane ]
     &          .or. ((molty(bi) .eq. 2) .and. (moltybj .eq. 6))        ! [ nonane-acetic_acid ]
     &          .or. ((molty(bi) .eq. 6) .and. (moltybj .eq. 2))) then  ! [ acetic_acid-nonane ]
	         enecut=ecwn/temp                                       ! $ecwn energy value is used here (recall this was set in fort.10).
	        elseif (((molty(bi) .eq. 2) .and. (moltybj .eq. 3))
     &           .or. ((molty(bi) .eq. 3) .and. (moltybj .eq. 2))
     &           .or. ((molty(bi) .eq. 2) .and. (moltybj .eq. 8))       ! [ nonane-octanol ]
     &           .or. ((molty(bi) .eq. 8) .and. (moltybj .eq. 2))) then ! [ octanol-nonane ]
	         enecut=ecnb/temp                                       ! $ecnb energy value is used here.
	        else                                                    
	         enecut=ecba/temp                                       ! $ecba energy value is used here.
	        endif                                                   
   	                                                                
	        ene=0d0                                                 
	        do runbi=1,bli                                          ! Each atom ($runbi) of molecule $bi compared with 
	         do runbj=1,blj                                         ! each atom ($runbj) of molecule $bj.
	          xij=x(bi,molty(bi),runbi)-x(ndx,moltybj,runbj)      
	          yij=y(bi,molty(bi),runbi)-y(ndx,moltybj,runbj)      
	          zij=z(bi,molty(bi),runbi)-z(ndx,moltybj,runbj)
	          disij2=xij**2+yij**2+zij**2
	          disij=disij2**0.5d0
	          epsij=(eps(molty(bi),runbi)*eps(moltybj,runbj))**0.5d0  ! Lorentz-Berthelot combining rules.
	          sigij=0.5d0*(sig(molty(bi),runbi)+sig(moltybj,runbj))   ! Lorentz-Berthelot combining rules.
	          sig6=(sigij/disij)**6
	          lenn=4*epsij*sig6*(sig6-1d0)                            ! Lennard-Jones potential.
	          coul1=qq(molty(bi),runbi)*qq(moltybj,runbj)/disij
	          coul=138.93545*coul1                                    ! Coulumbic potential.
	          ene=ene+lenn+coul                                       ! Total intermolecular interaction potential.
	         enddo
	        enddo                                     

   	        if ( ene .lt. enecut ) then          ! Alg found a match! Molecule passed criteria. 
   	         nuc(ndx)=frag                       ! Molecule $ndx assigned the unique cluster i.d.
	         vapor(ndx)=.FALSE.                  ! Molecule $ndx not a vapor as well since it has been clustered. 
	         pool(ndx)=.FALSE.                   ! FALSE, because molecule $ndx already part of a cluster.
	         RNDTAG(ndx)=.TRUE.                  ! TRUE, so later we know which molecules will star properties be applied to. 
	         clus(ndx)=.TRUE.                    ! Molecule $ndx tagged as satellite.  
	         mercnt=mercnt+1                     ! We add one to current mer count to register newly found molecule in cluster with star.
	         countb2=countb2+1                   ! Total number of satellites. 
   	         clus(bi)=.FALSE.                    ! FALSE, to remove satellite tag in current star section. 
	        endif            
	       endif            

        !  [ QUESTION ] clus(bi)? 
        !
        !  In a similar block in the star section, clus(bi) was not present. Is this a mistake?
        !  That's because that part was the star (never underwent clus(bi)=TRUE. 
        !  Star satellite should no longer be part of satellites as it has already been made a point of reference ("star")

  300          CONTINUE                              ! End point for Stillinger distance failed.
	      enddo                                  ! End scan of molecules in search box

	      if ( countb2 .eq. 0 ) then             ! Alg did not find a match. None passed Stillinger distance & energy criteria                          
	       clus(bi)=.FALSE.                      ! Satellite finished being a star, satellite tag can be turned off.
	      endif
	     endif                                  
	    enddo                                    ! End of branch analysis.

!-------III.F. DECIDING WHETHER TO FIND NEW STAR OR REITERATE BRANCHING

        !  [ PURPOSE ] In this section, 
        ! 
        !  Alg will evaluate whether it will repeat the searching procedure or continue to property analysis.
        !  This is an extremely significant section of the algorithm. If you think about it, 
        !  how would you have coded it, analysis by branching. If you're not careful you'd 
        !  end up thinking of a loop within a loop within a loop... 
        !  But alg only needed one extra loop. The key was the satellite i.d. 
        !  By using this tag, and removing the i.d. from satellites that have been analyzed,
        !  you can keep iterating back to the start (before analysis of
        !  first satellite) and keep using the satellite ids

	    iterate=.FALSE.   
	    do cc=1,pop                              ! After satellite is made a star, notice earlier that it's satellite tag is removed. 
	     if ( clus(cc) .eqv. .TRUE. ) then       ! This allows alg to determine if a satellite is yet to be analyzed by using the condition in this row.
	      iterate=.TRUE.                         ! If a satellite is found, iterate is set to TRUE.
	      EXIT
	     endif
	    enddo

	    if ( iterate .eqv. .TRUE. ) then         ! And if iterate is set to true, we return to start of branch analysis.
	     GOTO 100
	    else                                     ! If there are no more satellites left, we do this: 
	     do cc=1,pop                             ! Start tagging descriptors to molecules part of cluster.
	      if ( RNDTAG(cc) .eqv. .TRUE. ) then    ! RNDTAG will determine whether molecule is part of cluster.
	       mer(cc)=mercnt                        ! The total value counted by the monomer counter is associated to each molecule (cluster size)
	       RNDTAG(cc)=.FALSE.                    ! After which, RNDTAG is removed in preparation for next round.
	       tally(mercnt)=tally(mercnt)+1         ! Counts how many molecules have a size of $mercnt-mer
	      endif                                  
	     enddo
	     GOTO 200
	    endif

!-------III.G. TOP BLOCK TERMINATOR

  200       CONTINUE                                 ! End point for Stillinger energy failed. No branching employed.
	   endif                                     ! End point for empty tolerance. No branching employed.
 
  400      CONTINUE                                  ! End point for empty tolerance. No branching employed.
    	  endif                                      ! Molecule was in pool.
    	 enddo                                       ! All molecules in systems has been analyzed.
                                                     ! End of main detection scheme.
!-------IV. BEGIN CLUSTER RESTORATION

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#    PERFORM CLUSTER RESTORATION   #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"
          
!-------IV.A. COUNT RESTORED MOLECULES

        !  [ PURPOSE ] In this section, 
        ! 
        !  Alg will utilize the molecules restored earlier (start of Alg),
        !  as markers to determine likely places where a cluster was cut.
        !  The restored molecule is translated a box length in the axis where it's atom is protruding the system box.
        !  and from there, searching procedure will commence to test if molecules are close to it
        !  and whether it is enough to be clustered with it.

    	 counttransfer=0
   	 do pi=1,pop
	  pool(pi)=.TRUE.                            ! Returning all molecules back to pool. 
	  RNDTAG(pi)=.FALSE.                         ! Turning of all RNDTAGs.
    	  if( pbccheck(pi) .eqv. .TRUE. ) then       ! Checking for protruding molecules.
    	   counttransfer=counttransfer+1             
    	   pct(counttransfer)=pi                     ! All molecule i.d. of protruding molecules are associated with an additional i.d. 
    	  endif         
 	 enddo
c	 write(*,*) "Connecting clusters, total ",counttransfer

!-------IV.B. TRANSLATE RESTORED MOLECULES
       
   	 recov=0
    	 empfg=0
    	 do api=1,counttransfer                      ! Iterating through protruding molecules.
    	  pi=pct(api)                                ! We use new i.d. to call molecule i.d.
	  pool(pi)=.FALSE.                          
   	
	  if ( molty(pi) .eq. 1 ) then               ! Check molecule identity
	   lp=4                                      ! Assign length of molecule
	  elseif ( molty(pi) .eq. 2 ) then
	   lp=9
	  elseif ( molty(pi) .eq. 3 ) then
	   lp=6
	  elseif ( molty(pi) .eq. 4 ) then
	   lp=5
    	  elseif ( molty(pi) .eq. 5 ) then
   	   lp=3
    	  elseif ( molty(pi) .eq. 6 ) then
   	   lp=5
    	  elseif ( molty(pi) .eq. 7 ) then
   	   lp=1
    	  elseif ( molty(pi) .eq. 8 ) then
   	   lp=10
	  endif
    	               
    	  if ( pxmore(pi) .eqv. .TRUE. ) then        ! Endif is not used because protrusion may be in more than one axis.
    	   sux=0d0
    	   do pb=1,lp
    	    x(pi,molty(pi),pb)=x(pi,molty(pi),pb)-sybx
    	    sux=sux+x(pi,molty(pi),pb)
   	   enddo
    	  endif
    	  if ( pymore(pi) .eqv. .TRUE. ) then
    	   suy=0d0
    	   do pb=1,lp
    	    y(pi,molty(pi),pb)=y(pi,molty(pi),pb)-sybx
    	    suy=suy+y(pi,molty(pi),pb)
   	   enddo
    	  endif
   	  if ( pzmore(pi) .eqv. .TRUE. ) then
    	   suz=0d0
    	   do pb=1,lp
    	    z(pi,molty(pi),pb)=z(pi,molty(pi),pb)-sybx
    	    suz=suz+z(pi,molty(pi),pb)
    	   enddo 
    	  endif
    	  if ( pxless(pi) .eqv. .TRUE. ) then
    	   sux=0d0
   	   do pb=1,lp
   	    x(pi,molty(pi),pb)=x(pi,molty(pi),pb)+sybx
    	    sux=sux+x(pi,molty(pi),pb)
    	   enddo
    	  endif
    	  if ( pyless(pi) .eqv. .TRUE. ) then
    	   suy=0d0
   	   do pb=1,lp
    	    y(pi,molty(pi),pb)=y(pi,molty(pi),pb)+sybx
    	    suy=suy+y(pi,molty(pi),pb)
    	   enddo
    	  endif
    	  if ( pzless(pi) .eqv. .TRUE. ) then
    	   suz=0d0
    	   do pb=1,lp
    	    z(pi,molty(pi),pb)=z(pi,molty(pi),pb)+sybx  
    	    suz=suz+z(pi,molty(pi),pb)
   	   enddo  
    	  endif

!-------IV.C. REGISTER CENTER OF GRAVITY 

	  cgx(pi)=sux/DBLE(l) 
	  cgy(pi)=suy/DBLE(l) 
	  cgz(pi)=suz/DBLE(l)

!-------V. BEGIN STILLINGER TEST FOR RESTORATION 

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#    STILLINGER TEST RESTORATION   #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"

!-------V.A. GENERATE SEARCH BOX AROUND STAR MOLECULE

        ! [ PURPOSE ] In this section, 
        ! 
        ! searching procedure is commenced again to search for molecules within vicinity
        ! and test whether clustering is present.

	  RNDTAG(pi)=.TRUE.            
	  mnx=cgx(pi)-tolp             ! Only molecules within tolerance will be evaluated against Stillinger distance & energy criteria.
	  mny=cgy(pi)-tolp
	  mnz=cgz(pi)-tolp
	  mxx=cgx(pi)+tolp
	  mxy=cgy(pi)+tolp
	  mxz=cgz(pi)+tolp

!-------V.B. SCAN SEARCH BOX

	  countp1=0                                                   ! Counter for number of molecules within tolerance. 
	  do pj=1,pop                                                 ! Compare translated molecule to molecule pool.  
	   if ( (nuc(pi) .ne. nuc(pj))                                ! Some molecules that were translated may have been part of same cluster,       
     &     .and. ((cgx(pj) .gt. mnx) .and. (cgx(pj) .lt. mxx))        ! only after checking do we check if molecule $pj within tolerance.
     &     .and. ((cgy(pj) .gt. mny) .and. (cgy(pj) .lt. mxy))
     &     .and. ((cgz(pj) .gt. mnz) .and. (cgz(pj) .lt. mxz)) ) then 
	    countp1=countp1+1                                
	    pndxx(countp1)=pj                                         ! Associating a new i.d. for translated molecule to original molecule i.d.
	   endif                                                     
	  enddo                                                      

	  if (countp1 .eq. 0 ) then                                   ! Alg did not find a match. No molecule found within tolerance.
	   GOTO 401                                                   ! No need to evaluate Stillinger distance & energy.  
	  else                                                       
                                                                     
!-------V.C. LOOP ALL MOLECULES IN POOL3 (1 to countp1)               
                                                                     
	   do pnj=1,countp1                                           ! Iterate through all molecules inside tolerance. 
	    pndx=pndxx(pnj)                                            
	    moltypj=molty(pndx)                                      
                                                                     
!-------V.D. TEST NONBONDED PAIRWISE INTERACTIONS                    
                                                                     
	    if ( molty(pi) .eq. 1 ) then                              ! Calling the identity of molecule $pi
	     pli=4                                                    ! 1 =SOL; 2 =NON; 3=BUT; 4 =NH3 
	    elseif ( molty(pi) .eq. 2 ) then
	     pli=9
	    elseif ( molty(pi) .eq. 3 ) then
	     pli=6
	    elseif ( molty(pi) .eq. 4 ) then
	     pli=5
    	    elseif ( molty(pi) .eq. 5 ) then
    	     pli=3
    	    elseif ( molty(pi) .eq. 6 ) then
    	     pli=5
    	    elseif ( molty(pi) .eq. 7 ) then
    	     pli=1
    	    elseif ( molty(pi) .eq. 8 ) then
    	     pli=10
	    endif

	    if ( moltypj .eq. 1) then                                 ! Calling identity of molecule $pj.
	     plj=4
	    elseif ( moltypj .eq. 2 ) then
	     plj=9
	    elseif ( moltypj .eq. 3 ) then
	     plj=6
	    elseif ( moltypj .eq. 4 ) then
	     plj=5
    	    elseif ( moltypj .eq. 5 ) then
    	     plj=3
    	    elseif ( moltypj .eq. 6 ) then
    	     plj=5
    	    elseif ( moltypj .eq. 7 ) then
    	     plj=1
    	    elseif ( moltypj .eq. 8 ) then
    	     plj=10
	    endif

	    nullbox=.TRUE.                                   
	    do runpi=1,pli                                            ! Each atom ($runpi) of molecule $pi compared with 
	     do runpj=1,plj                                           ! each atom ($runpj) of molecule $pj.
	      sigij=0.5d0*(sig(molty(pi),runpi)+sig(moltypj,runpj))   ! Lorentz-Berthelot combining rules
	      discut=cutoff*sigij                                     ! Distance cut-off length
	      xij=x(pi,molty(pi),runpi)-x(pndx,moltypj,runpj)
	      yij=y(pi,molty(pi),runpi)-y(pndx,moltypj,runpj)
	      zij=z(pi,molty(pi),runpi)-z(pndx,moltypj,runpj)
	      disij2=xij**2+yij**2+zij**2
   	      disij=disij2**0.5d0
	      if ( disij .lt. discut ) then                           ! Checking whether molecule is within Stillinger distance cutoff.
   	       nullbox=.FALSE.                                        ! If FALSE, distance test passed.  
	       EXIT                                                   ! No need to test other atoms. Proceed to Stillinger energy criterion.
	      endif   
	     enddo                                          

	     if ( nullbox .eqv. .FALSE. ) then                         
	      EXIT                                      
	     endif
	    enddo                                     
	  
	    if ( nullbox .eqv. .TRUE. ) then                          ! Stillinger distance failed. Skip Stillinger energy test (essentially proceeding to next molecule in search box). 
	     GOTO 41                                                  ! Even without this line, a conditional is present before energy test that makes sure only 
	    endif                                                     ! those that passed may enter energy test.
                                                             
	    if ( disij .lt. discut ) then                             ! Commence stillinger energy test for molecules that passed distance criterion.
	      if (((molty(pi) .eq. 1) .and. (moltypj .eq. 2))         ! We need to call molecule identity because energy criteria will be different for each 
     &       .or. ((molty(pi) .eq. 2) .and. (moltypj .eq. 1))         ! molecule reacting
     &       .or. ((molty(pi) .eq. 2) .and. (moltypj .eq. 4))
     &       .or. ((molty(pi) .eq. 4) .and. (moltypj .eq. 2))
     &       .or. ((molty(pi) .eq. 2) .and. (moltypj .eq. 5))         ! [ nonane-methanol ] 
     &       .or. ((molty(pi) .eq. 5) .and. (moltypj .eq. 2))         ! [ methanol-nonane ] 
     &       .or. ((molty(pi) .eq. 2) .and. (moltypj .eq. 6))         ! [ nonane-acetic_acid ]
     &       .or. ((molty(pi) .eq. 6) .and. (moltypj .eq. 2)))  then  ! [ acetic_acid-nonane ]
	       enecut=ecwn/temp                                       ! $ecwn energy value is used here (recall this was set in fort.10).
	     elseif (((molty(pi) .eq. 2) .and. (moltypj .eq. 3))
     &          .or. ((molty(pi) .eq. 3) .and. (moltypj .eq. 2))
     &          .or. ((molty(pi) .eq. 2) .and. (moltypj .eq. 8))
     &          .or. ((molty(pi) .eq. 8) .and. (moltypj .eq. 2))) then
    	       enecut=ecnb/temp                                       ! $ecba energy value is used here.
	     else
	       enecut=ecba/temp
	     endif
	
	     ene=0d0
	     do runpi=1,pli                                           ! Each atom ($runpi) of molecule $pi compared with 
	      do runpj=1,plj                                          ! each atom ($runpj) of molecule $pj.
	       xij=x(pi,molty(pi),runpi)-x(pndx,moltypj,runpj)        
	       yij=y(pi,molty(pi),runpi)-y(pndx,moltypj,runpj)        
	       zij=z(pi,molty(pi),runpi)-z(pndx,moltypj,runpj)
	       disij2=xij**2+yij**2+zij**2
	       disij=disij2**0.5d0
	       epsij=(eps(molty(pi),runpi)*eps(moltypj,runpj))**0.5d0 ! Lorentz-Berthelot combining rules
	       sigij=0.5d0*(sig(molty(pi),runpi)+sig(moltypj,runpj))  ! Lorentz-Berthelot combining rules
	       sig6=(sigij/disij)**6  
	       lenn=4*epsij*sig6*(sig6-1d0)                           ! Lennard-Jones potential
	       coul1=qq(molty(pi),runpi)*qq(moltypj,runpj)/disij      
	       coul=138.93545*coul1                                   ! Coulumbic potential
	       ene=ene+lenn+coul                                      ! Total intermolecular interaction potential
    	      enddo
    	     enddo    

        !  [ QUESTION ] Why no branching algorithm? 
        ! 
        !  This is another crucial part of the algorithm. Translated molecule only needs to match with one molecule.
        !  Then using the unique cluster i.d., it can look for all other molecules clustered to that one molecule it found.

!-------V.E. RECOVERING CUT CLUSTER

        !  [ PURPOSE ] In this section, 
        !          
        !  all molecules connected to $pi and $pndx will have their descriptors be the same.
        !  whichever has the lower cluster i.d. between cluster A (where
        !  molecule $pi is part of) and cluster B (where molecule $pndx
        !  is part of) will be the point of reference. Size of each cluster is added and all molecules 
        !  are linked to that new cluster size. Tally of molecules that were initially recorded for 
        !  ${mercnt}-mer are also deducted accordingly.

	     if (ene .lt. enecut) then             ! Alg found match! Molecule passed all criteria.                   
    	      empfg=empfg+1                        
    	      recov=recov+1                                               
    	      fragp=nuc(pi)                        ! calling unique cluster i.d. of $pi ($pi is current star translated molecule) 
    	      fragpp=nuc(pndx)                     ! calling unique cluster i.d. of $pndx ($pndx is current molecule in search box that passed Stillinger distance and energy criteria)
    	      do ppp=1,pop                         
    	       fragppp=nuc(ppp)                    ! calling unique cluster i.d. of all other molecules
    	       if ( (fragppp .eq. fragp) 
     &         .or. (fragppp .eq. fragpp) )then 
    	        RNDTAG(ppp)=.TRUE.                 ! recall that RNDTAG is needed whenever we will set descriptors. 
    	       endif                             
    	       if ( fragppp .eq. fragpp ) then     ! removing from pool all molecules in same cluster with molecule $pndx.     
    	        pool(ppp)=.FALSE.
    	       endif
    	      enddo
             
    	      merp=mer(pndx)+mer(pi)               ! Combining mercnt of $pi & $pndx.
   	      do prr=1,pop                         ! Replacing nuc to lower frag.
    	       if ( RNDTAG(prr) .eqv. .TRUE. ) then
    	        tally(mer(prr))=tally(mer(prr))-1 
    	        mer(prr)=merp
   	        tally(merp)=tally(merp)+1
    	        RNDTAG(prr)=.FALSE.
    	        if ( nuc(pi) .lt. nuc(pndx) ) then
    	         nuc(prr)=nuc(pi)
    	        else
    	         nuc(prr)=nuc(pndx)
    	        endif
    	       endif 
   	      enddo 
  	      EXIT
   	       
   	     else                         
!  	     write(*,*) "MERGING FAILED",INT(ene*100/enecut),"%"
	     endif
	    endif                      

  41        CONTINUE    ! End point for Stillinger distance failed. Try another molecule in tolerance.
    	   enddo        ! End scan of molecules in tolerance.                
    	  endif         

  401     CONTINUE                                 ! End point for empty tolerance.
   	  call UNDO_TRANSLATE(pi)                  ! to avoid any complications, the translated molecules are returned
   	 enddo                                     ! end of iterating through all protruding molecules                         

        !  [ OPTIMIZATION CHECK ] Why undo the translation? 
        ! 
        !  This was our line of thinking. More than two molecules from the same cluster could be protruding. Let us label them molecule 124, and 723. 
        !  Molecule 124 is translated first, searches it's vicinity and undergoes restoration. Next, 723 is translated. Our worry was that 723 might see 124.
        !  But thinking about it now, it shouldn't have been a problem. Part of the conditional when matching is that cluster i.d. wasn't the same. 
        !  There shouldn't have been any problems right?  

!-------V.F. REORGANIZING FRAGS 

	 if ( recov .ne. 0 ) then                  ! Reorganizing frag
! 	  write(*,*)"RECOV",recov                  ! Removing gaps within frag count
 	  do frags=3,frag                          ! All frags are initially tagged as zero. 
 	   nonzero(frags)=.FALSE.
 	  enddo
                                                   ! Purpose of this block is to look for a non-zero benchmark.
   	 do fragb=3,frag                           ! Once the benchmark is found, we mark the frag right before it.
   	  do ppi=1,pop                             ! 1 is our lowest unique cluster i.d., so obviously it can never be a benchmark. 2 cannot be a benchmark because 1 can never be 0.  
   	   if ( nuc(ppi) .eq. fragb ) then         ! So we start with 3 in (fragb=3,frag).       
   	    gap=0                                  ! We can determine a non-zero by evaluating all molecules against a unique cluster i.d.
   	    nonzero(fragb)=.TRUE.                  ! If no molecules are associated with it, it is a zero frag.
    	    fragbb=fragb-1                         
    	    EXIT                                   
   	   endif
   	  enddo

    	  if ( nonzero(fragb) .eqv. .TRUE. ) then  ! It is important to make sure that current fragb (our benchmark frag) is non-zero.                
   	   do ppf=2,fragbb                         
   	    if ( nonzero(ppf) .eqv. .FALSE. ) then ! Searching for zero (nonzero=FALSE) (remember all i.d.s were initially set with a zero tag).
   	     nogap=.FALSE.                         ! Initially setting nogap as FALSE (same as saying gap = true).
   	     do ppj=1,pop                          ! Go through all molecules. 
   	      if ( nuc(ppj) .eq. ppf ) then        ! See if i.d. of any molecule is equal to current star $ppf.
   	       nogap=.TRUE.                        ! If an i.d. was found, $ppf is not a zero frag, so we make nogap=.TRUE. 
   	       nonzero(ppf)=.TRUE.                 ! and change it's nonzero status to TRUE.
   	      endif
   	     enddo                                 ! After going through all molecules, and the nogap=FALSE and nonzero=FALSE status never changed,
                                                   ! we are sure that a gap is present, that ppf is a zero frag.
	     if ( nogap .eqv. .FALSE. ) then       
	      gap=gap+1                            ! If the condition was met, we add 1 to gap counter.
	     endif
	    endif
	   enddo                                   ! We move on to next $ppf.  
	   if ( gap .ne. 0 ) then                  ! Check if gap not equal to 0 (if we found some zero frags).
    	    newfrag=fragb-gap                      ! If yes, we subtract the value of benchmark i.d. by the gap.
    	    nonzero(fragb)=.FALSE.                 ! Interchange tagging of old frag (benchmark frag before subtracting) and new frag.   
    	    nonzero(newfrag)=.TRUE.
   	    do ppk=1,pop
    	     if ( nuc(ppk) .eq. fragb ) then       ! Finally, interchange cluster i.d.s.
  	      nuc(ppk)=newfrag 
    	     endif
   	    enddo
    	   endif
  	  endif 
   	 enddo
    	 frag=frag-empfg                           ! Subtracting total number of frags.
	 endif 

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#   CLUSTER DETECTION COMPLETED    #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#      jgcabinta.rtabag.2022       #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################" 

!-------VI. BEGIN PROPERTY AND BEHAVIOR ANALYSIS 

!	    write(*,*) "####################################"
!	    write(*,*) "#                                  #"
!	    write(*,*) "#          ANALYSIS REPORT         #"
!	    write(*,*) "#                                  #"
!	    write(*,*) "####################################"

!-------VI.A.1. SEARCHING FOR CUT CLUSTERS 

        ! PURPOSE: In this section, 
        ! 
        ! we restore clusters. Recall, the very first restoration was of cut molecules.
        ! The second restoration was "restoring clusters". However it didn't follow the same 
        ! logic as the first. In the first one, all atoms were translated. In the second, just one 
        ! molecule was translated. By manipulating the cluster i.d.s, the cluster was restored
        ! by attributes only. If it were to be visualized, the cluster was still cut. This was 
        ! sufficient as we only needed to correct cluster i.d., size, and size tally. After this, 
        ! the translated molecule was returned to avoid double reading. Now, for the third restoration,
        ! we need to translate all molecules. Many of the analysis schemes will require measuring the 
        ! diameter of the sphere, hence a true restoration is needed. 

 	 do ffi=1,frag      
 	  REVD(ffi)=.FALSE. 
 	  do ppi=1,pop
 	   if (( nuc(ppi) .eq. ffi ) .and. ( mer(ppi) .eq. 1 )) then    ! No need to check monomer frags.
 	    REVD(ffi)=.TRUE.                                            ! Tagged as reviewed if monomer found.
 	    EXIT
 	   endif
  	  enddo
  	 enddo
         
  	 do ffi=1,frag                                                  ! Going through all frags.
   	  if ( REVD(ffi) .eqv. .FALSE. ) then                           ! Checking for frags not yet reviewed.
   	   REVD(ffi)=.TRUE.
 	   do ppi=1,pop                                                 ! Choosing a reference molecule.
   	    if ( nuc(ppi) .eq. ffi ) then
    	     EXIT
    	    endif
  	   enddo
                                                                        ! Comparing all other molecules to chosen reference molecule.
 	   do ppj=1,pop                                                 ! Specifically, only molecules clustered to reference. 
    	    if ((nuc(ppj) .eq. ffi) .and. (ppj .ne. ppi)) then          ! Second condition makes sure that reference molecule is not compared to itself. 
	     if ( molty(ppj) .eq. 1 ) then                              ! Calling identity of accepted molecule $ppj.
 	      l=4
	     elseif ( molty(ppj) .eq. 2 ) then
 	      l=9
	     elseif ( molty(ppj) .eq. 3 ) then
 	      l=6
	     elseif ( molty(ppj) .eq. 4 ) then
 	      l=5
	     elseif ( molty(ppj) .eq. 5 ) then
 	      l=3
	     elseif ( molty(ppj) .eq. 6 ) then
 	      l=5
	     elseif ( molty(ppj) .eq. 8 ) then
 	      l=10
  	     endif
  	     difx=x(ppi,molty(ppi),1)-x(ppj,molty(ppj),1)               ! Calculating distance between first atom of reference molecule, 
   	     dify=y(ppi,molty(ppi),1)-y(ppj,molty(ppj),1)               ! and first atom of accepted molecule.
 	     difz=z(ppi,molty(ppi),1)-z(ppj,molty(ppj),1)

!-------VI.A.2. RECOVERING CUT CLUSTER

    	     if ( abs(difx) .gt. sybx1 ) then                           ! Translate if clusters deemed to be cut.
    	      if ( x(ppj,molty(ppj),1) .gt. sybx1 ) then                ! Depending on whether it is greater or less than sybx1,
    	       do h=1,l                                                 ! translate by a system box length all the atoms of 
 	        x(ppj,molty(ppj),h)=x(ppj,molty(ppj),h)-sybx            ! accepted molecule.
    	       enddo
    	      else
    	       do h=1,l
    	        x(ppj,molty(ppj),h)=x(ppj,molty(ppj),h)+sybx
    	       enddo
    	      endif
    	     endif
   	
    	     if ( abs(dify) .gt. sybx1 ) then                          ! For y-axis.
    	      if ( y(ppj,molty(ppj),1) .gt. sybx1 ) then
    	       do h=1,l
   	        y(ppj,molty(ppj),h)=y(ppj,molty(ppj),h)-sybx
    	       enddo
    	      else
 	       do h=1,l
 	        y(ppj,molty(ppj),h)=y(ppj,molty(ppj),h)+sybx
    	       enddo
    	      endif
    	     endif
  	
    	     if ( abs(difz) .gt. sybx1 ) then                          ! For z-axis.
    	      if ( z(ppj,molty(ppj),1) .gt. sybx1 ) then
    	       do h=1,l
    	        z(ppj,molty(ppj),h)=z(ppj,molty(ppj),h)-sybx
    	       enddo
    	      else
    	       do h=1,l
    	        z(ppj,molty(ppj),h)=z(ppj,molty(ppj),h)+sybx
    	       enddo
    	      endif
    	     endif
    	
    	    endif
    	   enddo
  	  endif
 	 enddo

!-------VI.B.1. RADIAL DISTRIBUTION PROFILE #1 - MATCHING TARGET TIME

!  |  1: WHOLE SYSTEM / A: MOLECULES vs. RADIUS  |
!  |  1: WHOLE SYSTEM / B: DENSITY   vs. RADIUS  |

        ! [ PURPOSE ]
        !
        ! The scope of this profile is the entire system. r=0 is set at the  
        ! geometric center of all the molecules in the system. 
        ! Bins always begin at the center.
        ! If increment set to 3Å
        !   bin_1: r= 0Å to r= 3Å ; V=(4/3)π(3)^3  ; V=  113.097 Å^3 
        !   bin_2: r= 0Å to r= 6Å ; V=(4/3)π(6)^3  ; V=  904.779 A^3
        !   bin_3: r= 0Å to r= 9Å ; V=(4/3)π(9)^3  ; V= 3053.623 Å^3
        !
        ! The number of molecules in each bin is counted, and 
        ! divided by the volume of each bin to obtain the density.  
        !   bin_1: n=  5 ; D=  5/ 113.097 Å^3 = 0.0442098 Å^-3
        !   bin_2: n= 12 ; D= 12/ 904.779 A^3 = 0.0132629 A^-3
        !   bin_3: n= 18 ; D= 18/3053.623 Å^3 = 0.0058946 Å^-3

 	 if ( RPRD1 .eq. 1 ) then           
 	 qik = ctme - (1d-5)
 	 qok = ctme + (1d-5)
 	 
 	 do np = 1,5                        ! roll call of target time.
 	  if ( ( qik .lt. gt( np ) )  .and.
     &         ( qok .gt. gt( np ) ) ) then

!-------VI.B.2. RADIAL DISTRIBUTION PROFILE #1 - PREPARING REPORT FILES

 	   dp1h1 = '(a1,a13,f10.5)'
 	   dp1h2 = '(a1,a14,a17,a10,a19,2a17,a25,a15)' 
 	   dp1h3 = '(a1,a8,a17,a22,8(a23,a18))' 
 	   dp1h4 = '(a1,a5,3(a8,f8.3))'
 	   dp1c1 = '(i9,i18,i13,f21.3,2f17.3,f23.3,f15.3)'
 	   dp1c2 = '(i7,f17.3,f22.3,8(i23,f18.8))'
 	   
 	   qlv = 2000
 	   vvv = qlv + np
 	
 	   rmd = 3000
 	   rmc = rmd + np
 	
 	   lrb = 350             ! total r bins.
 	   bnr = 0.3d0           ! bin interval (3 Å).
 	
 	   write( vvv,dp1h1 ) "#","TIME (ns) =",ctme
 	   write( vvv,dp1h1 ) "#"
 	   write( rmc,dp1h1 ) "#","TIME (ns) =",ctme
 	   write( rmc,dp1h1 ) "#"
 	
 	   write( rmc,dp1h3 ) "#","BIN ID","RADIUS (A)","VOLUME (A^3)",
     &                       "N [ALL]", "DENSITY (A^-3)",
     &                       "N [WAT]", "DENSITY (A^-3)",
     &                       "N [NON]", "DENSITY (A^-3)",
     &                       "N [BUT]", "DENSITY (A^-3)",
     &                       "N [AMM]", "DENSITY (A^-3)",
     &                       "N [MET]", "DENSITY (A^-3)",
     &                       "N [ACT]", "DENSITY (A^-3)",
     &                       "N [OCT]", "DENSITY (A^-3)"

!-------VI.B.3. RADIAL DISTRIBUTION PROFILE #1 - SETTING STORAGE VARIABLES

 	   do rrf = 1,lrb                   ! collects value in chosen time step.
 	    genmol( np,rrf ) = 0            ! counts all molecules in bin r.
 	    sumwat( np,rrf ) = 0            ! only counts water in bin r.
 	    sumnon( np,rrf ) = 0            ! only counts nonane in bin r.
 	    sumbut( np,rrf ) = 0            ! only counts butanol in bin r.
 	    sumamm( np,rrf ) = 0            ! only counts ammonia in bin r.
 	    summet( np,rrf ) = 0            ! only counts methanol in bin r.
 	    sumact( np,rrf ) = 0            ! only counts acetate in bin r.
 	    sumoct( np,rrf ) = 0            ! only counts octanol in bin r.
 	   enddo

!-------VI.B.4. RADIAL DISTRIBUTION PROFILE #1 - SEARCHING CENTER r=0

  	   do ppi = 1,pop
  	    typcnt1 = 0             
  	    typcnt2 = 0             
  	    typcnt3 = 0             
  	    typcnt4 = 0                
  	    typcnt5 = 0                
  	    typcnt6 = 0                
  	    typcnt8 = 0                
  	    coox    = 0d0                    ! coordinates for center.
  	    cooy    = 0d0
  	    cooz    = 0d0
  	    if     ( molty( ppi ) .eq. 1 ) then
  	     typcnt1 = typcnt1 + 1
  	     l = 4
  	    elseif ( molty( ppi ) .eq. 2 ) then
  	     typcnt2 = typcnt2 + 1
  	     l = 9
  	    elseif ( molty( ppi ) .eq. 3 ) then
  	     typcnt3 = typcnt3 + 1
  	     l = 6
  	    elseif ( molty( ppi ) .eq. 4 ) then
  	     typcnt4 = typcnt4 + 1
  	     l = 5
  	    elseif ( molty( ppi ) .eq. 5 ) then      
  	     typcnt5 = typcnt5 + 1
  	     l = 3
  	    elseif ( molty( ppi ) .eq. 6 ) then      
  	     typcnt6 = typcnt6 + 1
  	     l = 5
  	    elseif ( molty( ppi ) .eq. 8 ) then      
  	     typcnt8 = typcnt8 + 1
  	     l = 10
  	    endif
  	    do k = 1,l                                              
  	     coox = coox + x( ppi,molty( ppi ),k )   
  	     cooy = cooy + y( ppi,molty( ppi ),k )    
  	     cooz = cooz + z( ppi,molty( ppi ),k )   
  	    enddo                                                
  	   enddo
  	
  	   coll=DBLE(typcnt1*4+typcnt2*9+
     &               typcnt3*6+typcnt4*5+
     &               typcnt5*3+typcnt6*5+
     &               typcnt8*10)
  	   cddx=coox/coll
  	   cddy=cooy/coll
  	   cddz=cooz/coll
  	
 	   write( vvv,dp1h4 ) "#","R=0","X:",cddx,"Y:",cddy,"Z:",cddz
 	   write( vvv,dp1h2 ) "#"
 	   write( vvv,dp1h2 ) "#","MOLECULE NO.","MOLECULE TYPE","ATOM",
     &                        "X COORD","Y COORD","Z COORD",
     &                        "DIST. F.C.","REP. DIST."    ! DISTANCE FROM CENTER, REPRESENTATIVE DISTANCE.

!-------VI.B.5. RADIAL DISTRIBUTION PROFILE #1 - SEARCHING MOLEC DISTANCE FROM CENTER

   	   do ppk = 1,pop                                       
    	    cox = 0d0                                           
    	    coy = 0d0                                           
    	    coz = 0d0   
 	    if     ( molty( ppk ) .eq. 1 ) then
 	     li = 4
 	    elseif ( molty( ppk ) .eq. 2 ) then
 	     li = 9
 	    elseif ( molty( ppk ) .eq. 3 ) then
 	     li = 6
 	    elseif ( molty( ppk ) .eq. 4 ) then
 	     li = 5
 	    elseif ( molty( ppk ) .eq. 5 ) then
 	     li = 3
 	    elseif ( molty( ppk ) .eq. 6 ) then
 	     li = 5
 	    elseif ( molty( ppk ) .eq. 8 ) then
 	     li = 10
 	    endif
 	
 	    write(vvv,*) " "
 	
 	    do kk = 1,li
 	     cox = x( ppk,molty( ppk ),kk )  
 	     coy = y( ppk,molty( ppk ),kk )    
 	     coz = z( ppk,molty( ppk ),kk )
	     rsx  = cddx - cox                                  
	     rsy  = cddy - coy                                          
	     rsz  = cddz - coz                                          
 	     rad2 = SQRT( rsx**2 + rsy**2 + rsz**2 )
 	
 	     if ( kk .eq. 1 ) then               ! point of reference of molec is atom closest to center.
 	      brad2 = rad2 
 	     else
 	      if ( rad2 .lt. brad2 ) then
 	       brad2 = rad2
 	      endif
 	     endif
 	
 	  write(vvv,dp1c1) ppk,molty(ppk),kk,cox,coy,coz,rad2,brad2
 	    enddo

!-------VI.B.6. RADIAL DISTRIBUTION PROFILE #1 - STORING DIST BET CELL CENTER AND INDIV MOLEC

	    rrb  = 0d0
	    do rrf = 1,lrb                       ! records dist of molec from cell center.                              
	     rrb = rrb + bnr                                     
	     if ( brad2 .le. rrb ) then 
	      genmol( np,rrf ) = genmol( np,rrf ) + 1                  
 	      if    (  molty( ppk ) .eq. 1 ) then                     
 	       sumwat( np,rrf ) = sumwat( np,rrf ) + 1
 	       EXIT                              ! once recorded, proceeds to next molec.                             
 	      elseif ( molty( ppk ) .eq. 2 ) then                   
 	       sumnon( np,rrf ) = sumnon( np,rrf ) + 1           
 	       EXIT                                               
 	      elseif ( molty( ppk ) .eq. 3 ) then                     
 	       sumbut( np,rrf ) = sumbut( np,rrf ) + 1           
 	       EXIT                                               
 	      elseif ( molty( ppk ) .eq. 4 ) then                   
 	       sumamm( np,rrf)  = sumamm( np,rrf ) + 1           
 	       EXIT                                               
 	      elseif ( molty( ppk ) .eq. 5 ) then                   
 	       summet( np,rrf)  = summet( np,rrf ) + 1           
 	       EXIT                                               
 	      elseif ( molty( ppk ) .eq. 6 ) then                   
 	       sumact( np,rrf)  = sumact( np,rrf ) + 1           
 	       EXIT                                               
 	      elseif ( molty( ppk ) .eq. 8 ) then                   
 	       sumoct( np,rrf)  = sumoct( np,rrf ) + 1           
 	       EXIT                                               
 	      endif
 	     endif                               ! only molec within radius bin may procee.
 	    enddo                                ! end of iterating radius bins.
 	   enddo                                 ! end of iterating all molec in current time step.
   	   
   	   do rrf=1,lrb                          ! including inner bin r in each bin.
   	    rrg=rrf-1
   	    genmol( np,rrf ) = genmol( np,rrf ) + genmol( np,rrg )
   	    sumwat( np,rrf ) = sumwat( np,rrf ) + sumwat( np,rrg )
   	    sumnon( np,rrf ) = sumnon( np,rrf ) + sumnon( np,rrg )
   	    sumbut( np,rrf ) = sumbut( np,rrf ) + sumbut( np,rrg )
   	    sumamm( np,rrf ) = sumamm( np,rrf ) + sumamm( np,rrg )
   	    summet( np,rrf ) = summet( np,rrf ) + summet( np,rrg )
   	    sumact( np,rrf ) = sumact( np,rrf ) + sumact( np,rrg )
   	    sumoct( np,rrf ) = sumoct( np,rrf ) + sumoct( np,rrg )
   	   enddo
   
!-------VI.B.7. RADIAL DISTRIBUTION FILE #1 - REMOVING EXTRA VALUES

   	   cgen=genmol(np,lrb)                   ! remove excess bins (bins are repeated at the end).
   	   do rrf=lrb,1,-1                       ! reverse increment until current value is different from last value.
   	    if (genmol(np,rrf) .ne. cgen) then   
   	     rrh=rrf+1
 	     psk=rrh
   	     genmol(np,rrh)=kksg                 ! the true last value is also deleted, needs this to return.
   	     sumwat(np,rrh)=kksw
   	     sumnon(np,rrh)=kksn
   	     sumbut(np,rrh)=kksb
   	     sumamm(np,rrh)=kksa
   	     summet(np,rrh)=kksm
   	     sumact(np,rrh)=kksc
   	     sumoct(np,rrh)=kkso
   	     EXIT
   	    else
   	     kksg=genmol(np,rrf)
   	     kksw=sumwat(np,rrf)
   	     kksn=sumnon(np,rrf)
   	     kksb=sumbut(np,rrf)
   	     kksa=sumamm(np,rrf)
   	     kksm=summet(np,rrf)
   	     kksc=sumact(np,rrf)
   	     kkso=sumoct(np,rrf)
   	     genmol(np,rrf)=0d0                  ! set excess bins to 0.
   	     sumwat(np,rrf)=0d0
   	     sumnon(np,rrf)=0d0
   	     sumbut(np,rrf)=0d0
   	     sumamm(np,rrf)=0d0
   	     summet(np,rrf)=0d0
   	     sumact(np,rrf)=0d0
   	     sumoct(np,rrf)=0d0
   	    endif 
   	   enddo 

!-------VI.B.8. RADIAL DISTRIBUTION PROFILE #1 - CALCULATING VOLUME AND DENSITY

 	   rrb=0d0
 	   do rrf = 1,psk
 	    rrb = rrb + bnr  
 	    divvol = (4d0/3d0)*3.1415927d0*(rrb*10)**3         ! volume of bin r.
 	
 	    ndgen( np,rrf ) = DBLE(genmol( np,rrf )) / divvol              
 	    ndwat( np,rrf ) = DBLE(sumwat( np,rrf )) / divvol               
 	    ndnon( np,rrf ) = DBLE(sumnon( np,rrf )) / divvol 
 	    ndbut( np,rrf ) = DBLE(sumbut( np,rrf )) / divvol 
 	    ndamm( np,rrf ) = DBLE(sumamm( np,rrf )) / divvol 
 	    ndmet( np,rrf ) = DBLE(summet( np,rrf )) / divvol 
 	    ndact( np,rrf ) = DBLE(sumact( np,rrf )) / divvol 
 	    ndoct( np,rrf ) = DBLE(sumoct( np,rrf )) / divvol 
 	    
 	    der = rrb - bnr
 	    write(rmc,dp1c2) rrf,der*10d0,divvol,
     &          genmol( np,rrf ),ndgen( np,rrf ),
     &          sumwat( np,rrf ),ndwat( np,rrf ),
     &          sumnon( np,rrf ),ndnon( np,rrf ),
     &          sumbut( np,rrf ),ndbut( np,rrf ),
     &          sumamm( np,rrf ),ndamm( np,rrf ),
     &          summet( np,rrf ),ndmet( np,rrf ),
     &          sumact( np,rrf ),ndact( np,rrf ),
     &          sumoct( np,rrf ),ndoct( np,rrf )
   	   enddo
  	  endif                              ! only time matching with target time may proceed.
  	 enddo                               ! roll call of target time.
  	 endif
 
!-------VI.C.1. RADIAL DISTRIBUTION PROFILE #2 - MATCHING TARGET TIME AND MER

!  |  2: CLUSTER OR CLUSTER RANGE  / A: MOLECULES vs. RADIUS  |
!  |  2: CLUSTER OR CLUSTER RANGE  / B: DENSITY   vs. RADIUS  |

        ! [ PURPOSE ] 
        !
        ! The scope of this profile is a single cluster. r=0 is set at the 
        ! geometric center of cluster. Bins always begin at center.
        ! If increment set to 3Å
        !   bin_1: r= 0Å to r= 3Å ; V=(4/3)π(3)^3  ; V=  113.097Å^3 
        !   bin_2: r= 0Å to r= 6Å ; V=(4/3)π(6)^3  ; V=  904.779A^3
        !   bin_3: r= 0Å to r= 9Å ; V=(4/3)π(9)^3  ; V= 3053.623Å^3
        !
        ! The number of molecules in each bin is counted, and 
        ! divided by the volume of each bin to obtain the density.  
        !   bin_1: n=  5 ; D=  5/ 113.097 Å^3 = 0.0442098 Å^-3
        !   bin_2: n= 12 ; D= 12/ 904.779 A^3 = 0.0132629 A^-3
        !   bin_3: n= 18 ; D= 18/3053.623 Å^3 = 0.0058946 Å^-3

 	 if ( RPRD2 .eq. 1 ) then 
 	 qik = ctme - (1d-5)
 	 qok = ctme + (1d-5)

 	 do np = 1,5                         ! roll call of target time.
 	  if ( ( qik .lt. gt( np ) )  .and.
     &         ( qok .gt. gt( np ) ) ) then
 	   do nth = 1,le( np )               ! roll call of target mers.
 	    cclm(nth)  = 0                   ! to count all mers present in range.
 	    crm( nth ) = 0                   ! counts total frag in chosen range.
 	    nin = gmin( np,nth )             ! lower limit of chosen range.
 	    nfn = gmax( np,nth )             ! upper limit of chosen range.
 	    nap = snap( np,nth )             ! target mer in range.

!-------VI.C.2. RADIAL DISTRIBUTION PROFILE #2 - PREPARING REPORT FILES

 	    dp1h1 = '(a1,a13,f10.5)'
 	    dp1h2 = '(a1,a14,a17,a10,a19,2a17,a25,a15)' 
 	    dp1h3 = '(a1,a8,a17,a22,8(a23,a18))' 
 	    dp1h4 = '(a1,a5,3(a8,f8.3))' 
 	    dp1h5 = '(a1,a17,i8,a8,i8)'
 	    dp1c1 = '(i9,i18,i13,f21.3,2f17.3,f23.3,f15.3)'
 	    dp1c2 = '(i7,f17.3,f22.3,8(i23,f18.8))'
 	    dp2h3 = '(a1,a11,i8,a14,i8)' 
 	    
 	    qlv = 2006
 	    mxv = np*3
 	    vvv = qlv + mxv + nth
 	    
 	    rmd = 3006
 	    mxv = np*3
 	    rmc = rmd + mxv + nth
 	    
 	    lrb=350             !total r bins
 	    bnr=0.3d0           !bin interval
 	    
 	    write( vvv,dp1h1 ) "#","TIME (ns) =",ctme
 	    write( vvv,dp1h1 ) "#"
 	    write( rmc,dp1h1 ) "#","TIME (ns) =",ctme
 	    write( rmc,dp1h1 ) "#"
 	    write( rmc,dp1h5 ) "#","MER RANGE =",nin,"to",nfn
 	    write( rmc,dp1h5 ) "#"        
 	
    	    do ffi = 1,frag
    	     FRGPOOL( ffi ) = .TRUE.
    	    enddo

!-------VI.C.3. RADIAL DISTRIBUTION PROFILE #2 - SETTING STORAGE VARIABLES

 	    do rrf = 1,lrb                  ! collects value for cluster sizes in chosen range.
 	     genmol2( nth,rrf ) = 0
 	     sumwat2( nth,rrf ) = 0      
 	     sumnon2( nth,rrf ) = 0
 	     sumbut2( nth,rrf ) = 0
 	     sumamm2( nth,rrf ) = 0
 	     summet2( nth,rrf ) = 0
             sumact2( nth,rrf ) = 0
 	     sumoct2( nth,rrf ) = 0
 	    enddo
 
!-------VI.C.4. RADIAL DISTRIBUTION PROFILE #2 - ITERATE ALL MER SIZE
 
   	    do mmi  = 1,pop           
 	     if ( ( tally( mmi ) .ne. 0   ) .and. 
     &            ( mmi          .ge. nin ) .and. 
     &            ( mmi          .le. nfn ) ) then
  	      tfrag = tally( mmi )/mmi
   	            
 	      do rrf = 1,lrb                  ! collects value for clusters of same size.   
 	       rmmg2( rrf ) = 0d0
 	       rmmw2( rrf ) = 0d0         
 	       rmmn2( rrf ) = 0d0
 	       rmmb2( rrf ) = 0d0
 	       rmma2( rrf ) = 0d0
 	       rmmm2( rrf ) = 0d0
 	       rmmc2( rrf ) = 0d0
 	       rmmo2( rrf ) = 0d0
 	      enddo
 
!-------VI.C.5. RADIAL DISTRIBUTION PROFILE #2 - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE
 
    	      do ppi = 1,pop  
 	       if ( mer( ppi ) .eq. mmi ) then
 	        do ffj = 1,frag  
 	         if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &                ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	      
 	          FRGPOOL( ffj ) = .FALSE.  
 	          typcnt1 = 0             
 	          typcnt2 = 0             
 	          typcnt3 = 0             
 	          typcnt4 = 0                
 	          typcnt5 = 0                
 	          typcnt6 = 0                
 	          typcnt8 = 0                
 	          coox    = 0d0              ! coordinates for cluster center.
 	          cooy    = 0d0
 	          cooz    = 0d0
 	      
 	          do rrf = 1,lrb             ! collects value for present cluster.     
 	           rgen2( rrf ) = 0
 	           rwat2( rrf ) = 0
 	           rnon2( rrf ) = 0
 	           rbut2( rrf ) = 0
 	           ramm2( rrf ) = 0
 	           rmet2( rrf ) = 0
 	           ract2( rrf ) = 0
 	           roct2( rrf ) = 0
 	          enddo
  	
 	   write( rmc,dp2h3 ) "#","MER =",mmi,"FRAG",ffj
 	   write( rmc,dp2h3 ) "#"
                          
!-------VI.C.6. RADIAL DISTRIBUTION PROFILE #2 - INDEXING MOLECULES INVOLVED IN LOCATED FRAG
 
 	          do ppj = 1,pop                                               
 	           if ( nuc( ppj ) .eq. ffj ) then
 	            if     ( molty( ppj ) .eq. 1 ) then
 	             typcnt1 = typcnt1 + 1
 	             l = 4
 	            elseif ( molty( ppj ) .eq. 2 ) then
 	             typcnt2 = typcnt2 + 1
 	             l = 9
 	            elseif ( molty( ppj ) .eq. 3 ) then
 	             typcnt3 = typcnt3 + 1
 	             l = 6
 	            elseif ( molty( ppj ) .eq. 4 ) then
 	             typcnt4 = typcnt4 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 5 ) then
 	             typcnt5 = typcnt5 + 1
 	             l = 3
 	            elseif ( molty( ppj ) .eq. 6 ) then
 	             typcnt6 = typcnt6 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 8 ) then
 	             typcnt8 = typcnt8 + 1
 	             l = 10
 	            endif
 
!-------VI.C.7. RADIAL DISTRIBUTION PROFILE #2 - SEARCHING CENTER OF WHOLE CLUSTER
 
 	            do k = 1,l                                              
 	             coox = coox + x( ppj,molty( ppj ),k )   
 	             cooy = cooy + y( ppj,molty( ppj ),k )    
 	             cooz = cooz + z( ppj,molty( ppj ),k )   
  	            enddo                                                
 	           endif                                                 
 	          enddo
 	
 	          coll=DBLE(typcnt1*4+typcnt2*9+
     &                      typcnt3*6+typcnt4*5+
     &                      typcnt5*3+typcnt6*5+
     &                      typcnt8*10)
  	          cddx=coox/coll
  	          cddy=cooy/coll
  	          cddz=cooz/coll
  	
 	   write( vvv,dp1h4 ) "#","R=0","X:",cddx,"Y:",cddy,"Z:",cddz
 	   write( vvv,dp1h2 ) "#"
 	   write( vvv,dp1h2 ) "#","MOLECULE NO.","MOLECULE TYPE","ATOM",
     &                        "X COORD","Y COORD","Z COORD",
     &                        "DIST. F.C.","REP. DIST."    ! DISTANCE FROM CENTER, REPRESENTATIVE DISTANCE.
 
!-------VI.C.8. RADIAL DISTRIBUTION PROFILE #2 - SEARCHING REFERENCE ATOM OF INDIV MOLEC
 
   	          do ppk = 1,pop                                       
  	           if ( nuc( ppk ) .eq. ffj ) then
  	            if     ( molty( ppk ) .eq. 1 ) then
    	             li = 4
 	            elseif ( molty( ppk ) .eq. 2 ) then
 	             li = 9
  	            elseif ( molty( ppk ) .eq. 3 ) then
 	             li = 6
 	            elseif ( molty( ppk ) .eq. 4 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 5 ) then
  	             li = 3
 	            elseif ( molty( ppk ) .eq. 6 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 8 ) then
  	             li = 10
 	            endif
 	          
 	            write(vvv,*) " "
 	
  	            do kk = 1,li
   	             cox = x( ppk,molty( ppk ),kk )   
   	             coy = y( ppk,molty( ppk ),kk )    
   	             coz = z( ppk,molty( ppk ),kk )   
 	             rsx  = cddx - cox                                  
 	             rsy  = cddy - coy                                          
 	             rsz  = cddz - coz                                          
 	             rad2 = SQRT( rsx**2 + rsy**2 + rsz**2 )
 	           
 	             if ( kk .eq. 1 ) then   ! point of reference of molec is atom closest to center.
 	              brad2 = rad2
 	             else
 	              if ( rad2 .lt. brad2 ) then
 	               brad2 = rad2
 	              endif
 	             endif
  	                    
 	      write(vvv,dp1c1) ppk,molty(ppk),kk,cox,coy,coz,rad2,brad2
 	            enddo
 
!-------VI.C.9. RADIAL DISTRIBUTION PROFILE #2 - DIST BET CLUSTER CENTER AND INDIV MOLEC
 
 	            rrb = 0d0
 	            do rrf = 1,lrb           ! records dist of molec from center.
 	             rrb = rrb + bnr                                     
 	             if ( brad2 .le. rrb ) then 
 	              rgen2( rrf ) = rgen2( rrf ) + 1                         
 	              if    (  molty( ppk ) .eq. 1 ) then                     
 	               rwat2( rrf ) = rwat2( rrf ) + 1
 	               EXIT                  ! once recorded, proceeds to next molec.                             
 	              elseif ( molty( ppk ) .eq. 2 ) then                   
 	               rnon2( rrf ) = rnon2( rrf ) + 1           
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 3 ) then                     
 	               rbut2( rrf ) = rbut2( rrf ) + 1           
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 4 ) then                   
 	               ramm2( rrf ) = ramm2( rrf ) + 1           
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 5 ) then                   
 	               rmet2( rrf ) = rmet2( rrf ) + 1 
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 6 ) then                   
 	               ract2( rrf ) = ract2( rrf ) + 1           
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 8 ) then                   
 	               roct2( rrf ) = roct2( rrf ) + 1           
 	               EXIT                                               
 	              endif
 	             endif                   ! only molec within radius bin may proceed.                             
 	            enddo                    ! end of iterating radius bins.
 	           endif                     ! only molec of current frag/cluster may proceed.                         
 	          enddo                      ! end of iterating all molec in current cluster.                                                   
 	
 	          do rrf = 1,lrb 
 	           rrg=rrf-1     
 	           rgen2( rrf ) = rgen2( rrf ) + rgen2( rrg ) 
 	           rwat2( rrf ) = rwat2( rrf ) + rwat2( rrg )      
 	           rnon2( rrf ) = rnon2( rrf ) + rnon2( rrg )      
 	           rbut2( rrf ) = rbut2( rrf ) + rbut2( rrg )      
 	           ramm2( rrf ) = ramm2( rrf ) + ramm2( rrg )      
 	           rmet2( rrf ) = rmet2( rrf ) + rmet2( rrg )      
 	           ract2( rrf ) = ract2( rrf ) + ract2( rrg )      
 	           roct2( rrf ) = roct2( rrf ) + roct2( rrg )      
 	          enddo 
                  
!-------VI.C.10. RADIAL DISTRIBUTION PROFILE #2 - SUMMING COUNT FOR CLUSTERS OF SAME SIZE

 	          do rrf = 1,lrb
 	           rmmg2( rrf ) = rmmg2( rrf ) + DBLE(rgen2( rrf )) 
 	           rmmw2( rrf ) = rmmw2( rrf ) + DBLE(rwat2( rrf ))      
 	           rmmn2( rrf ) = rmmn2( rrf ) + DBLE(rnon2( rrf ))      
 	           rmmb2( rrf ) = rmmb2( rrf ) + DBLE(rbut2( rrf ))      
 	           rmma2( rrf ) = rmma2( rrf ) + DBLE(ramm2( rrf ))
 	           rmmm2( rrf ) = rmmm2( rrf ) + DBLE(rmet2( rrf ))
 	           rmmc2( rrf ) = rmmc2( rrf ) + DBLE(ract2( rrf ))
 	           rmmo2( rrf ) = rmmo2( rrf ) + DBLE(roct2( rrf ))
 	          enddo                                                    
 	         endif                       ! only molec contained in current frag/cluster may proceed.
 	        enddo                        ! end of searching frags with current cluster size.
 	       endif                         ! only molec contained in current mer size may proceed.
 	      enddo                          ! end of iterating through all molecs.
 
!-------VI.C.11. RADIAL DISTRIBUTION PROFILE #2 - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
 
 	      do rrf = 1,lrb
 	       rmmg2( rrf ) = rmmg2( rrf ) / DBLE( tfrag )
 	       rmmw2( rrf ) = rmmw2( rrf ) / DBLE( tfrag )
 	       rmmn2( rrf ) = rmmn2( rrf ) / DBLE( tfrag )
 	       rmmb2( rrf ) = rmmb2( rrf ) / DBLE( tfrag )
 	       rmma2( rrf ) = rmma2( rrf ) / DBLE( tfrag )
 	       rmmm2( rrf ) = rmmm2( rrf ) / DBLE( tfrag )
 	       rmmc2( rrf ) = rmmc2( rrf ) / DBLE( tfrag )
 	       rmmo2( rrf ) = rmmo2( rrf ) / DBLE( tfrag )
 	      enddo
 
!-------VI.C.12. RADIAL DISTRIBUTION PROFILE #2 - SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
 
 	      crm( nth ) = crm( nth ) + 1
    	      do rrf = 1,lrb
 	       genmol2( nth,rrf ) = genmol2( nth,rrf ) + rmmg2( rrf ) 
    	       sumwat2( nth,rrf ) = sumwat2( nth,rrf ) + rmmw2( rrf )     
    	       sumnon2( nth,rrf ) = sumnon2( nth,rrf ) + rmmn2( rrf )
    	       sumbut2( nth,rrf ) = sumbut2( nth,rrf ) + rmmb2( rrf )
    	       sumamm2( nth,rrf ) = sumamm2( nth,rrf ) + rmma2( rrf )
    	       summet2( nth,rrf ) = summet2( nth,rrf ) + rmmm2( rrf )
    	       sumact2( nth,rrf ) = sumact2( nth,rrf ) + rmmc2( rrf )
    	       sumoct2( nth,rrf ) = sumoct2( nth,rrf ) + rmmo2( rrf )
    	      enddo
  	     endif
    	    enddo                            !end of iterating all cluster sizes

 	   write( rmc,dp1h3 ) "#","BIN ID","RADIUS (A)","VOLUME (A^3)",
     &                       "N [ALL]", "DENSITY (A^-3)",
     &                       "N [WAT]", "DENSITY (A^-3)",
     &                       "N [NON]", "DENSITY (A^-3)",
     &                       "N [BUT]", "DENSITY (A^-3)",
     &                       "N [AMM]", "DENSITY (A^-3)",
     &                       "N [MET]", "DENSITY (A^-3)",
     &                       "N [ACT]", "DENSITY (A^-3)",
     &                       "N [OCT]", "DENSITY (A^-3)"
            
!-------VI.C.13. RADIAL DISTRIBUTION PROFILE #2 - REMOVING EXTRA VALUES

 	    cgen=genmol2(nth,lrb)          
 	    do rrf=lrb,1,-1                       ! reverse increment until current value is different from last value.
   	     if ( ( genmol2( nth,rrf ) .ne. cgen ) .or.
     &            ( rrf .eq. 1) ) then
   	      rrh=rrf+1     
   	
 	      if ( cclm(nth) .eq. 0 ) then     
 	       cclm(nth) = rrh                    ! cclm counts all mers present in range
 	      elseif ( rrh .gt. cclm(nth) ) then
 	       cclm(nth) = rrh
 	      endif
 	
 	      if (rrf .ne. 1) then
   	       genmol2( nth,rrh )=kksg            ! the true last value is also deleted, needs this to return.
   	       sumwat2( nth,rrh )=kksw
   	       sumnon2( nth,rrh )=kksn
   	       sumbut2( nth,rrh )=kksb
   	       sumamm2( nth,rrh )=kksa
   	       summet2( nth,rrh )=kksm
   	       sumact2( nth,rrh )=kksc
   	       sumoct2( nth,rrh )=kksc
   	       EXIT
   	      else                                ! this fixes the issue if only one bin is occupied.
   	       cclm(nth)=cclm(nth)-1
   	       EXIT
   	      endif
   	
   	     else
   	      kksg=genmol2( nth,rrf )
   	      kksw=sumwat2( nth,rrf )
   	      kksn=sumnon2( nth,rrf )
   	      kksb=sumbut2( nth,rrf )
   	      kksa=sumamm2( nth,rrf )
   	      kksm=summet2( nth,rrf )
   	      kksc=sumact2( nth,rrf )
   	      kkso=sumoct2( nth,rrf )
   	      genmol2( nth,rrf )=0d0              ! set excess bins to 0.
   	      sumwat2( nth,rrf )=0d0
   	      sumnon2( nth,rrf )=0d0
   	      sumbut2( nth,rrf )=0d0
   	      sumamm2( nth,rrf )=0d0
   	      summet2( nth,rrf )=0d0
   	      sumact2( nth,rrf )=0d0
   	      sumoct2( nth,rrf )=0d0
   	     endif 
   	    enddo 

!-------VI.C.14. RADIAL DISTRIBUTION PROFILE #2 - AVERAGING COUNT TO REPRESENT CHOSEN RANGE
            
 	    rrb=0d0
 	    do rrf = 1,cclm(nth)
 	    rrb = rrb + bnr  
 	     if ( crm( nth ) .eq. 0 ) then
 	      dvgen( nth,rrf ) = 0d0
 	      dvwat( nth,rrf ) = 0d0
 	      dvnon( nth,rrf ) = 0d0
 	      dvbut( nth,rrf ) = 0d0
 	      dvamm( nth,rrf ) = 0d0
 	      dvmet( nth,rrf ) = 0d0
 	      dvact( nth,rrf ) = 0d0
 	      dvoct( nth,rrf ) = 0d0
 	     else
 	      dvgen( nth,rrf ) = genmol2(nth,rrf) /DBLE(crm(nth)) 
 	      dvwat( nth,rrf ) = sumwat2(nth,rrf) /DBLE(crm(nth))
 	      dvnon( nth,rrf ) = sumnon2(nth,rrf) /DBLE(crm(nth))
 	      dvbut( nth,rrf ) = sumbut2(nth,rrf) /DBLE(crm(nth))
 	      dvamm( nth,rrf ) = sumamm2(nth,rrf) /DBLE(crm(nth))
 	      dvmet( nth,rrf ) = summet2(nth,rrf) /DBLE(crm(nth))
 	      dvact( nth,rrf ) = sumact2(nth,rrf) /DBLE(crm(nth))
 	      dvoct( nth,rrf ) = sumoct2(nth,rrf) /DBLE(crm(nth))
 	     endif
 
!-------VI.C.15. RADIAL DISTRIBUTION PROFILE #2 - CALCULATING VOLUME AND DENSITY
 
 	     divvol = (4d0/3d0)*3.1415927d0*(rrb*10)**3         ! volume of bin r.
 	     ndgen( nth,rrf ) = dvgen( nth,rrf ) / divvol
 	     ndwat( nth,rrf ) = dvwat( nth,rrf ) / divvol !distribution function               
 	     ndnon( nth,rrf ) = dvnon( nth,rrf ) / divvol !distribution function 
 	     ndbut( nth,rrf ) = dvbut( nth,rrf ) / divvol !distribution function 
 	     ndamm( nth,rrf ) = dvamm( nth,rrf ) / divvol !distribution function 
 	     ndmet( nth,rrf ) = dvmet( nth,rrf ) / divvol !distribution function 
 	     ndact( nth,rrf ) = dvact( nth,rrf ) / divvol !distribution function 
 	     ndoct( nth,rrf ) = dvoct( nth,rrf ) / divvol !distribution function 
 	
 	     rmd = 3006
 	     mxv = np*3
 	     rmc = rmd + mxv + nth
 	     der = rrb - bnr
 	   
 	     write(rmc,dp1c2) rrf,der*10d0,divvol,
     &                 genmol2( nth,rrf ),ndgen( nth,rrf ),
     &                 sumwat2( nth,rrf ),ndwat( nth,rrf ),
     &                 sumnon2( nth,rrf ),ndnon( nth,rrf ),
     &                 sumbut2( nth,rrf ),ndbut( nth,rrf ),
     &                 sumamm2( nth,rrf ),ndamm( nth,rrf ), 
     &                 summet2( nth,rrf ),ndmet( nth,rrf ), 
     &                 sumact2( nth,rrf ),ndact( nth,rrf ), 
     &                 sumoct2( nth,rrf ),ndoct( nth,rrf ) 
   	    enddo
  	   enddo                             ! roll call of target cluster sizes.
  	  endif                              ! only time matching with target time may proceed. 
   	 enddo                               ! roll call of target time.
   	 endif

!-------VI.D.1. RADIAL DISTRIBUTION PROFILE #3 - MATCHING TARGET TIME AND MER

!  |  3: CLUSTER OR CLUSTER RANGE  / A: MOLECULES vs. RADIUS  |
!  |  3: CLUSTER OR CLUSTER RANGE  / B: DENSITY   vs. RADIUS  |

        ! [ PURPOSE ] 
        !
        ! The scope of this profile is a single cluster. r=0 is set at the 
        ! geometric center of cluster. Bins don't always begin from center. 
        ! If increment set to 3Å:
        !   bin_1: r= 0Å to r< 3Å ; V=(4/3)π ( (3)^3-(0)^3 )  ; V=  113.097Å^3 
        !   bin_2: r= 3Å to r< 6Å ; V=(4/3)π ( (6)^3-(3)^3 )  ; V=  791.681A^3
        !   bin_3: r= 6Å to r< 9Å ; V=(4/3)π ( (9)^3-(6)^3 )  ; V= 2148.849Å^3
        !
        ! The number of molecules in each bin is counted, and 
        ! divided by the volume of each bin to obtain the density.  
        !   bin_1: n=  5 ; D= 5/ 113.097 Å^3 = 0.0442098 Å^-3
        !   bin_2: n=  7 ; D= 7/ 791.681 A^3 = 0.0084194 A^-3
        !   bin_3: n=  6 ; D= 6/2148.849 Å^3 = 0.0027922 Å^-3

 	 if ( RPRD3 .eq. 1 ) then 
 	 qik = ctme - (1d-5)
 	 qok = ctme + (1d-5)
  	            
 	 do np = 1,5                         ! roll call of target time.
 	  if ( ( qik .lt. gt( np ) )  .and.
     &         ( qok .gt. gt( np ) ) ) then
 	   do nth = 1,le( np )               ! roll call of target mers
 	    cclm2( nth ) = 0 
 	    kin( nth ) = 0
 	    crm3( nth ) = 0                  ! counts total frag in chosen range.
 	    nin = gmin( np,nth )             ! lower limit of chosen range.
 	    nfn = gmax( np,nth )             ! upper limit of chosen range.
 	    nap = snap( np,nth )             ! target mer in range.
            
!-------VI.D.2. RADIAL DISTRIBUTION PROFILE #3 - PREPARING REPORT FILES

 	    dfh1 = '(a1,a17,f10.5)'
 	    dfh2 = '(a1,a17,i8,a8,i8)' 
 	    dfh3 = '(a1,a11,i8,a14,i8)' 
 	    dfh4 = '(a1,a17,a29,a12,a27,2a26)'
 	    dfh5 = '(a1,a13,f10.5)' 
 	    dfc1 = '(f21.16,f45.16,3f26.16)' 
 	    
  	    qlv = 2046
  	    mxv = np*3
  	    vvv = qlv + mxv + nth
  	
 	    rmd = 3046
 	    mxv = np*3
 	    rmc = rmd + mxv + nth
 	
 	    lrb=350             ! total r bins.
 	    bnr=0.3d0           ! bin interval.
 	
 	    write( vvv,dfh5 ) "#","TIME (ns) =",ctme
 	    write( vvv,dfh5 ) "#"
  	    write( rmc,dfh5 ) "#","TIME (ns) =",ctme
  	    write( rmc,dfh5 ) "#"
  	    write( rmc,dfh2 ) "#","MER RANGE =",nin,"to",nfn
  	    write( rmc,dfh2 ) "#"
  	
    	    do ffi = 1,frag
    	     FRGPOOL( ffi ) = .TRUE.
    	    enddo
 
!-------VI.D.3. RADIAL DISTRIBUTION PROFILE #3 - SETTING STORAGE VARIABLES

 	    do rrf = 1,lrb                    ! collects value for cluster sizes in chosen range.
 	     genmol3( nth,rrf ) = 0
 	     sumwat3( nth,rrf ) = 0      
 	     sumnon3( nth,rrf ) = 0
 	     sumbut3( nth,rrf ) = 0
 	     sumamm3( nth,rrf ) = 0
 	     summet3( nth,rrf ) = 0
 	     sumact3( nth,rrf ) = 0
 	     sumoct3( nth,rrf ) = 0
  	    enddo

!-------VI.D.4. RADIAL DISTRIBUTION PROFILE #3 - ITERATE ALL MER SIZE
 
   	    do mmi  = 1,pop           
 	     if ( ( tally( mmi ) .ne. 0   ) .and. 
     &            ( mmi          .ge. nin ) .and. 
     &            ( mmi          .le. nfn ) ) then
  	      tfrag = tally( mmi )/mmi
  	      lres=0d0  
  	          
 	      do rrf = 1,lrb                  ! collects value for clusters of same size.
  	       rmmg3( rrf ) = 0d0
 	       rmmw3( rrf ) = 0d0         
 	       rmmn3( rrf ) = 0d0
 	       rmmb3( rrf ) = 0d0
 	       rmma3( rrf ) = 0d0
 	       rmmm3( rrf ) = 0d0
 	       rmmc3( rrf ) = 0d0
 	       rmmo3( rrf ) = 0d0
 	      enddo
 
!-------VI.D.5. RADIAL DISTRIBUTION PROFILE #3 - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE
 
    	      do ppi = 1,pop  
 	       if ( mer( ppi ) .eq. mmi ) then
 	        do ffj = 1,frag  
 	         if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &                ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	      
 	          FRGPOOL( ffj ) = .FALSE.  
 	          typcnt1 = 0             
 	          typcnt2 = 0             
 	          typcnt3 = 0             
 	          typcnt4 = 0                
 	          typcnt5 = 0                
 	          typcnt6 = 0                
 	          typcnt8 = 0                
 	          coox    = 0d0              ! coordinates for cluster center.
 	          cooy    = 0d0
 	          cooz    = 0d0
 	      
 	          do rrf = 1,lrb             ! collects value for present cluster.     
 	           rgen3( rrf ) = 0
 	           rwat3( rrf ) = 0
 	           rnon3( rrf ) = 0
 	           rbut3( rrf ) = 0
 	           ramm3( rrf ) = 0
 	           rmet3( rrf ) = 0
 	           ract3( rrf ) = 0
 	           roct3( rrf ) = 0
 	          enddo
  	
 	          write( rmc,dfh3 ) "#","MER =",mmi,"FRAG",ffj
 	          write( rmc,dfh3 ) "#"
                          
!-------VI.D.6. RADIAL DISTRIBUTION PROFILE #3 - INDEXING MOLECULES INVOLVED IN LOCATED FRAGS
 
 	          do ppj = 1,pop                                               
 	           if ( nuc( ppj ) .eq. ffj ) then
 	            if     ( molty( ppj ) .eq. 1 ) then
 	             typcnt1 = typcnt1 + 1
 	             l = 4
 	            elseif ( molty( ppj ) .eq. 2 ) then
 	             typcnt2 = typcnt2 + 1
 	             l = 9
 	            elseif ( molty( ppj ) .eq. 3 ) then
 	             typcnt3 = typcnt3 + 1
 	             l = 6
 	            elseif ( molty( ppj ) .eq. 4 ) then
 	             typcnt4 = typcnt4 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 5 ) then
 	             typcnt5 = typcnt5 + 1
 	             l = 3
 	            elseif ( molty( ppj ) .eq. 6 ) then
 	             typcnt6 = typcnt6 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 8 ) then
 	             typcnt8 = typcnt8 + 1
 	             l = 10
 	            endif
 
!-------VI.D.7. RADIAL DISTRIBUTION PROFILE #3 - SEARCHING CENTER OF WHOLE CLUSTER
 
 	            do k = 1,l                                              
 	             coox = coox + x( ppj,molty( ppj ),k )   
 	             cooy = cooy + y( ppj,molty( ppj ),k )    
 	             cooz = cooz + z( ppj,molty( ppj ),k )   
  	            enddo                                                
 	           endif                                                 
 	          enddo
 	
 	          coll=DBLE(typcnt1*4+typcnt2*9+
     &                      typcnt3*6+typcnt4*5+
     &                      typcnt5*3+typcnt6*5+
     &                      typcnt8*10)
  	          cddx=coox/coll
  	          cddy=cooy/coll
  	          cddz=cooz/coll
   	            
 	   write( vvv,dp1h4 ) "#","R=0","X:",cddx,"Y:",cddy,"Z:",cddz
 	   write( vvv,dp1h2 ) "#"
 	   write( vvv,dp1h2 ) "#","MOLECULE NO.","MOLECULE TYPE","ATOM",
     &                        "X COORD","Y COORD","Z COORD",
     &                        "DIST. F.C.","REP. DIST."    ! DISTANCE FROM CENTER, REPRESENTATIVE DISTANCE.

!-------VI.D.8. RADIAL DISTRIBUTION PROFILE #3 - SEARCHING REFERENCE ATOM OF INDIV MOLEC
 
   	          do ppk = 1,pop                                       
    	           cox = 0d0                                           
    	           coy = 0d0                                           
    	           coz = 0d0   
    	                 
  	           if ( nuc( ppk ) .eq. ffj ) then
  	            if     ( molty( ppk ) .eq. 1 ) then
    	             li = 4
 	            elseif ( molty( ppk ) .eq. 2 ) then
 	             li = 9
  	            elseif ( molty( ppk ) .eq. 3 ) then
 	             li = 6
 	            elseif ( molty( ppk ) .eq. 4 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 5 ) then
  	             li = 3
 	            elseif ( molty( ppk ) .eq. 6 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 8 ) then
  	             li = 8
 	            endif
 	
 	            write(vvv,*) " "
 	          
 	            do kk = 1,li
 	             cox = x( ppk,molty( ppk ),kk )   
 	             coy = y( ppk,molty( ppk ),kk )    
 	             coz = z( ppk,molty( ppk ),kk )   
 	             rsx  = cddx - cox                                  
 	             rsy  = cddy - coy                                          
 	             rsz  = cddz - coz                                          
 	             rad2 = SQRT( rsx**2 + rsy**2 + rsz**2 )
 	
 	             if ( kk .eq. 1 ) then   ! point of reference of molec is atom closest to center.
 	              brad2 = rad2
 	             else
 	              if ( rad2 .lt. brad2 ) then
 	               brad2 = rad2
 	              endif
 	             endif
  	                    
 	    write(vvv,dp1c1) ppk,molty(ppk),kk,cox,coy,coz,rad2,brad2
 	            enddo

!-------VI.D.9. RADIAL DISTRIBUTION PROFILE #3 - DIST BET CLUSTER CENTER AND INDIV MOLEC
 
 	            rrb = 0d0
 	            larg= 0d0                    
 	            do rrf = 1,lrb            ! records dist of molec from center.         
 	             rrb = rrb + bnr                                     
 	             if ( brad2 .le. rrb ) then   
 	              rgen3( rrf ) = rgen3( rrf ) + 1
 	              if    (  molty( ppk ) .eq. 1 ) then                     
 	               rwat3( rrf ) = rwat3( rrf ) + 1
 	               larg=rrb            
 	               EXIT                  ! once recorded, proceeds to next molec.
 	              elseif ( molty( ppk ) .eq. 2 ) then                   
 	               rnon3( rrf ) = rnon3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 3 ) then                     
 	               rbut3( rrf ) = rbut3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 4 ) then                   
 	               ramm3( rrf ) = ramm3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 5 ) then                   
 	               rmet3( rrf ) = rmet3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 6 ) then                   
 	               ract3( rrf ) = ract3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 8 ) then                   
 	               roct3( rrf ) = roct3( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              endif                                               
 	             endif                   ! only molec within radius bin may proceed.                             
 	            enddo                    ! end of iterating radius bins.
 	 
 	            if ( larg .gt. lres ) then
 	                lres = larg
 	            endif
 	 
 	           endif                     ! only molec of current frag/cluster may proceed.                         
 	          enddo                      ! end of iterating all molec in current cluster.                                                     
 
!-------VI.D.10. RADIAL DISTRIBUTION PROFILE #3 - SUMMING COUNT FOR CLUSTERS OF SAME SIZE
 
 	          do rrf = 1,lrb                                             
 	           rmmg3( rrf ) = rmmg3( rrf ) + DBLE(rgen3( rrf ))      
 	           rmmw3( rrf ) = rmmw3( rrf ) + DBLE(rwat3( rrf ))      
 	           rmmn3( rrf ) = rmmn3( rrf ) + DBLE(rnon3( rrf ))      
 	           rmmb3( rrf ) = rmmb3( rrf ) + DBLE(rbut3( rrf ))      
 	           rmma3( rrf ) = rmma3( rrf ) + DBLE(ramm3( rrf ))      
 	           rmmm3( rrf ) = rmmm3( rrf ) + DBLE(rmet3( rrf ))      
 	           rmmc3( rrf ) = rmmc3( rrf ) + DBLE(ract3( rrf ))      
 	           rmmo3( rrf ) = rmmo3( rrf ) + DBLE(roct3( rrf ))      
 	          enddo                                                     
 	         endif                       ! only molec contained in current frag/cluster may proceed.
 	        enddo                        ! end of searching frags with current cluster size.
 	       endif                         ! only molec contained in current mer size may proceed.
 	      enddo                          ! end of iterating through all molecs.
 
!-------VI.D.11. RADIAL DISTRIBUTION PROFILE #3 - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
 
 	      do rrf = 1,lrb
 	       rmmg3( rrf ) = rmmg3( rrf ) / DBLE( tfrag )
 	       rmmw3( rrf ) = rmmw3( rrf ) / DBLE( tfrag )
 	       rmmn3( rrf ) = rmmn3( rrf ) / DBLE( tfrag )
 	       rmmb3( rrf ) = rmmb3( rrf ) / DBLE( tfrag )
 	       rmma3( rrf ) = rmma3( rrf ) / DBLE( tfrag )
 	       rmmm3( rrf ) = rmmm3( rrf ) / DBLE( tfrag )
 	       rmmc3( rrf ) = rmmc3( rrf ) / DBLE( tfrag )
 	       rmmo3( rrf ) = rmmo3( rrf ) / DBLE( tfrag )
 	      enddo
 
!-------VI.D.12. RADIAL DISTRIBUTION PROFILE #3 - SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
 
 	      crm3( nth ) = crm3( nth ) + 1
    	      do rrf = 1,lrb
    	       genmol3( nth,rrf ) = genmol3( nth,rrf ) + rmmg3( rrf ) 
    	       sumwat3( nth,rrf ) = sumwat3( nth,rrf ) + rmmw3( rrf )     
    	       sumnon3( nth,rrf ) = sumnon3( nth,rrf ) + rmmn3( rrf )
    	       sumbut3( nth,rrf ) = sumbut3( nth,rrf ) + rmmb3( rrf )
    	       sumamm3( nth,rrf ) = sumamm3( nth,rrf ) + rmma3( rrf )
    	       summet3( nth,rrf ) = summet3( nth,rrf ) + rmmm3( rrf )
    	       sumact3( nth,rrf ) = sumact3( nth,rrf ) + rmmc3( rrf )
    	       sumoct3( nth,rrf ) = sumoct3( nth,rrf ) + rmmo3( rrf )
    	      enddo
  	     endif
    	    enddo                            ! end of iterating all cluster sizes.
  	           
 	    write( rmc,dp1h3 ) "#","BIN ID","RADIUS (A)","VOLUME (A^3)",
     &                       "N [ALL]", "DENSITY (A^-3)",
     &                       "N [WAT]", "DENSITY (A^-3)",
     &                       "N [NON]", "DENSITY (A^-3)",
     &                       "N [BUT]", "DENSITY (A^-3)",
     &                       "N [AMM]", "DENSITY (A^-3)",
     &                       "N [MET]", "DENSITY (A^-3)",
     &                       "N [ACT]", "DENSITY (A^-3)",
     &                       "N [OCT]", "DENSITY (A^-3)"
 
!-------VI.D.13. RADIAL DISTRIBUTION PROFILE #3 - REMOVING EXTRA VALUES	

 	    cclm2(nth) = 0
 	    do rrf = lrb,1,-1
 	     if ( genmol3(nth,rrf) .eq. 0 ) then
 	      cclm2(nth) = cclm2(nth) + 1
 	     else
 	      kin(nth) = lrb - cclm2(nth)
 	      EXIT
 	     endif
 	    enddo

!-------VI.D.14. RADIAL DISTRIBUTION PROFILE #3 - AVERAGING COUNT TO REPRESENT CHOSEN RANGE
 
 	    rrb=0d0
 	    do rrf = 1,kin(nth)
 	     rrb = rrb + bnr  
 	     if ( crm3( nth ) .eq. 0 ) then
 	      dvgen( nth,rrf ) = 0d0
 	      dvwat( nth,rrf ) = 0d0
 	      dvnon( nth,rrf ) = 0d0
 	      dvbut( nth,rrf ) = 0d0
 	      dvamm( nth,rrf ) = 0d0
 	      dvmet( nth,rrf ) = 0d0
 	      dvact( nth,rrf ) = 0d0
 	      dvoct( nth,rrf ) = 0d0
 	     else
 	      dvgen( nth,rrf ) = genmol3(nth,rrf)/DBLE(crm3(nth)) 
 	      dvwat( nth,rrf ) = sumwat3(nth,rrf)/DBLE(crm3(nth))
 	      dvnon( nth,rrf ) = sumnon3(nth,rrf)/DBLE(crm3(nth))
 	      dvbut( nth,rrf ) = sumbut3(nth,rrf)/DBLE(crm3(nth))
 	      dvamm( nth,rrf ) = sumamm3(nth,rrf)/DBLE(crm3(nth))
 	      dvmet( nth,rrf ) = summet3(nth,rrf)/DBLE(crm3(nth))
 	      dvact( nth,rrf ) = sumact3(nth,rrf)/DBLE(crm3(nth))
 	      dvoct( nth,rrf ) = sumoct3(nth,rrf)/DBLE(crm3(nth))
 	     endif
 
!-------VI.D.15. RADIAL DISTRIBUTION PROFILE #3 - CALCULATING VOLUME AND DENSITY
 
 	    divvol = (4d0/3d0)*3.1415927d0*(rrb*10)**3        
     &      - (4d0/3d0)*3.1415927d0*( rrb*10 - bnr*10 )**3 ! volume of bin r.
  	            
 	     
 	     ndgen( nth,rrf ) = dvgen( nth,rrf ) / divvol
 	     ndwat( nth,rrf ) = dvwat( nth,rrf ) / divvol                
 	     ndnon( nth,rrf ) = dvnon( nth,rrf ) / divvol  
 	     ndbut( nth,rrf ) = dvbut( nth,rrf ) / divvol  
 	     ndamm( nth,rrf ) = dvamm( nth,rrf ) / divvol  
 	     ndmet( nth,rrf ) = dvmet( nth,rrf ) / divvol  
 	     ndact( nth,rrf ) = dvact( nth,rrf ) / divvol  
 	     ndoct( nth,rrf ) = dvoct( nth,rrf ) / divvol  
 	     
 	     rmd = 3046
 	     mxv = np*3
 	     rmc = rmd + mxv + nth
 	     der = rrb - bnr
 	    
 	     write(rmc,dp1c2) rrf,der*10d0,divvol,
     &                 genmol3( nth,rrf ),ndgen( nth,rrf ),
     &                 sumwat3( nth,rrf ),ndwat( nth,rrf ),
     &                 sumnon3( nth,rrf ),ndnon( nth,rrf ),
     &                 sumbut3( nth,rrf ),ndbut( nth,rrf ),
     &                 sumamm3( nth,rrf ),ndamm( nth,rrf ), 
     &                 summet3( nth,rrf ),ndmet( nth,rrf ), 
     &                 sumact3( nth,rrf ),ndact( nth,rrf ), 
     &                 sumoct3( nth,rrf ),ndoct( nth,rrf ) 
   	    enddo
   	   enddo
  	  endif                              ! only time matching with target time may proceed. 
   	 enddo                               ! roll call of target time.
   	 endif

!-------VI.E.1. RADIAL DISTRIBUTION FUNCTION - MATCHING TARGET TIME AND MER

!  |   CLUSTER OR CLUSTER RANGE  / A: MOLECULES vs. RADIUS  |
!  |   CLUSTER OR CLUSTER RANGE  / B: g(r)      vs. RADIUS  |

        ! [ PURPOSE ] 
        !
        ! The scope of this profile is a single cluster. r=0 is set at the 
        ! geometric center of cluster. 
        ! Bins don't always begin from center. 
        ! If increment set to 3Å:
        !   bin_1: r= 0Å to r< 3Å ; V=(4/3)π ( (3)^3-(0)^3 )  ; V=  113.097Å^3 
        !   bin_2: r= 3Å to r< 6Å ; V=(4/3)π ( (6)^3-(3)^3 )  ; V=  791.681A^3
        !   bin_3: r= 6Å to r< 9Å ; V=(4/3)π ( (9)^3-(6)^3 )  ; V= 2148.849Å^3
        !
        ! The number of molecules in each bin is counted, and 
        ! divided by the volume of each bin to obtain the density.  
        !   bin_1: n=  5 ; D= 5/ 113.097 Å^3 = 0.0442098 Å^-3
        !   bin_2: n=  7 ; D= 7/ 791.681 A^3 = 0.0084194 A^-3
        !   bin_3: n=  6 ; D= 6/2148.849 Å^3 = 0.0027922 Å^-3
        ! 
        ! Normalization is employed by dividing with overall density.
        ! Overall Density = Overall Molecule Count / Overall Volume 
        !                 = 18 / 3053.623 Å^3
        !                 = 0.00589464 Å^-3
        ! g(r)= 0.0442098 Å^-3 /0.00589464 Å^-3  = 7.5 
        ! g(r)= 0.0084194 A^-3 /0.00589464 Å^-3  = 1.428315
        ! g(r)= 0.0027922 Å^-3 /0.00589464 Å^-3  = 0.473685

 	 if ( RDF .eq. 1 ) then 
 	 qik = ctme - (1d-5)
 	 qok = ctme + (1d-5)
  	            
 	 do np = 1,5                         ! roll call of target time.
 	  if ( ( qik .lt. gt( np ) )  .and.
     &         ( qok .gt. gt( np ) ) ) then
 	   do nth = 1,le( np )               ! roll call of target mers.
 	    cclm3( nth ) = 0 
 	    kin2( nth ) = 0
 	    crm4( nth ) = 0                  ! counts total frag in chosen range.
 	    nin = gmin( np,nth )             ! lower limit of chosen range.
 	    nfn = gmax( np,nth )             ! upper limit of chosen range.
 	    nap = snap( np,nth )             ! target mer in range.
            
!-------VI.E.2. RADIAL DISTRIBUTION FUNCTION - PREPARING REPORT FILES

 	    dfh1 = '(a1,a17,f10.5)'
 	    dfh2 = '(a1,a17,i8,a8,i8)' 
 	    dfh3 = '(a1,a11,i8,a14,i8)' 
 	    dfh4 = '(a1,a17,a29,a12,a27,2a26)'
 	    dfh5 = '(a1,a13,f10.5)' 
 	    dfc1 = '(f21.16,f45.16,3f26.16)' 
 	    dp1h3 = '(a1,a8,a17,a22,8(a23,a14))' 
 	    
  	    qlv = 2066
  	    mxv = np*3
  	    vvv = qlv + mxv + nth
  	
 	    rmd = 3066
 	    mxv = np*3
 	    rmc = rmd + mxv + nth
 	
 	    lrb=350             ! total r bins.
 	    bnr=0.3d0           ! bin interval.
 	
 	    write( vvv,dfh5 ) "#","TIME (ns) =",ctme
 	    write( vvv,dfh5 ) "#"
  	    write( rmc,dfh5 ) "#","TIME (ns) =",ctme
 	    write( rmc,dfh5 ) "#"
  	    write( rmc,dfh2 ) "#","MER RANGE =",nin,"to",nfn
  	    write( rmc,dfh2 ) "#"

    	    do ffi = 1,frag
    	     FRGPOOL( ffi ) = .TRUE.
    	    enddo
 
!-------VI.E.3. RADIAL DISTRIBUTION FUNCTION - SETTING STORAGE VARIABLES

 	    do rrf = 1,lrb                    ! collects value for cluster sizes in chosen range.
 	     genmol4( nth,rrf ) = 0
 	     sumwat4( nth,rrf ) = 0      
 	     sumnon4( nth,rrf ) = 0
 	     sumbut4( nth,rrf ) = 0
 	     sumamm4( nth,rrf ) = 0
 	     summet4( nth,rrf ) = 0
 	     sumact4( nth,rrf ) = 0
 	     sumoct4( nth,rrf ) = 0
  	    enddo

!-------VI.E.4. RADIAL DISTRIBUTION FUNCTION - ITERATE ALL MER SIZE
 
   	    do mmi  = 1,pop           
 	     if ( ( tally( mmi ) .ne. 0   ) .and. 
     &            ( mmi          .ge. nin ) .and. 
     &            ( mmi          .le. nfn ) ) then
  	      tfrag = tally( mmi )/mmi
  	      lres=0d0  
  	          
 	      do rrf = 1,lrb                  ! collects value for clusters of same size.
 	       rmmg4( rrf ) = 0d0
 	       rmmw4( rrf ) = 0d0         
 	       rmmn4( rrf ) = 0d0
 	       rmmb4( rrf ) = 0d0
 	       rmma4( rrf ) = 0d0
 	       rmmm4( rrf ) = 0d0
 	       rmmc4( rrf ) = 0d0
 	       rmmo4( rrf ) = 0d0
 	      enddo
 
!-------VI.E.5. RADIAL DISTRIBUTION FUNCTION - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE
 
    	      do ppi = 1,pop  
 	       if ( mer( ppi ) .eq. mmi ) then
 	        do ffj = 1,frag  
 	         if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &                ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	      
 	          FRGPOOL( ffj ) = .FALSE.  
 	          typcnt1 = 0             
 	          typcnt2 = 0             
 	          typcnt3 = 0             
 	          typcnt4 = 0                
 	          typcnt5 = 0                
 	          typcnt6 = 0                
 	          typcnt8 = 0                
 	          coox    = 0d0              ! coordinates for cluster center.
 	          cooy    = 0d0
 	          cooz    = 0d0
 	      
 	          do rrf = 1,lrb             ! collects value for present cluster.
 	           rgen4( rrf ) = 0
 	           rwat4( rrf ) = 0
 	           rnon4( rrf ) = 0
 	           rbut4( rrf ) = 0
 	           ramm4( rrf ) = 0
 	           rmet4( rrf ) = 0
 	           ract4( rrf ) = 0
 	           roct4( rrf ) = 0
 	          enddo
  	
 	          write( rmc,dfh3 ) "#","MER =",mmi,"FRAG",ffj
 	          write( rmc,dfh3 ) "#"
                          
!-------VI.E.6. RADIAL DISTRIBUTION FUNCTION - INDEXING MOLECULES INVOLVED IN LOCATED FRAG
 
 	          do ppj = 1,pop                                               
 	           if ( nuc( ppj ) .eq. ffj ) then
 	            if     ( molty( ppj ) .eq. 1 ) then
 	             typcnt1 = typcnt1 + 1
 	             l = 4
 	            elseif ( molty( ppj ) .eq. 2 ) then
 	             typcnt2 = typcnt2 + 1
 	             l = 9
 	            elseif ( molty( ppj ) .eq. 3 ) then
 	             typcnt3 = typcnt3 + 1
 	             l = 6
 	            elseif ( molty( ppj ) .eq. 4 ) then
 	             typcnt4 = typcnt4 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 5 ) then
 	             typcnt5 = typcnt5 + 1
 	             l = 3
 	            elseif ( molty( ppj ) .eq. 6 ) then
 	             typcnt6 = typcnt6 + 1
 	             l = 5
 	            elseif ( molty( ppj ) .eq. 8 ) then
 	             typcnt8 = typcnt8 + 1
 	             l = 5
 	            endif
 
!-------VI.E.7. RADIAL DISTRIBUTION FUNCTION - SEARCHING CENTER OF WHOLE CLUSTER
 
 	            do k = 1,l                                              
 	             coox = coox + x( ppj,molty( ppj ),k )   
 	             cooy = cooy + y( ppj,molty( ppj ),k )    
 	             cooz = cooz + z( ppj,molty( ppj ),k )   
  	            enddo                                                
 	           endif                                                 
 	          enddo
 	
 	          coll=DBLE(typcnt1*4+typcnt2*9+
     &                      typcnt3*6+typcnt4*5+
     &                      typcnt5*3+typcnt6*5+
     &                      typcnt8*10)
  	          cddx=coox/coll
  	          cddy=cooy/coll
  	          cddz=cooz/coll
   	            
 	   write( vvv,dp1h4 ) "#","R=0","X:",cddx,"Y:",cddy,"Z:",cddz
 	   write( vvv,dp1h2 ) "#"
 	   write( vvv,dp1h2 ) "#","MOLECULE NO.","MOLECULE TYPE","ATOM",
     &                        "X COORD","Y COORD","Z COORD",
     &                        "DIST. F.C.","REP. DIST."    ! DISTANCE FROM CENTER, REPRESENTATIVE DISTANCE.

!-------VI.E.8. RADIAL DISTRIBUTION FUNCTION - SEARCHING CENTER OF INDIV MOLEC
 
   	          do ppk = 1,pop                                       
    	           cox = 0d0                                           
    	           coy = 0d0                                           
    	           coz = 0d0   
    	                 
  	           if ( nuc( ppk ) .eq. ffj ) then
  	            if     ( molty( ppk ) .eq. 1 ) then
    	             li = 4
 	            elseif ( molty( ppk ) .eq. 2 ) then
 	             li = 9
  	            elseif ( molty( ppk ) .eq. 3 ) then
 	             li = 6
 	            elseif ( molty( ppk ) .eq. 4 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 5 ) then
  	             li = 3
 	            elseif ( molty( ppk ) .eq. 6 ) then
  	             li = 5
 	            elseif ( molty( ppk ) .eq. 8 ) then
  	             li = 10
 	            endif
 	
 	            write(vvv,*) " "
 	          
  	            do kk = 1,li
   	             cox = x( ppk,molty( ppk ),kk )   
   	             coy = y( ppk,molty( ppk ),kk )    
   	             coz = z( ppk,molty( ppk ),kk )   
 	             rsx  = cddx - cox                                  
 	             rsy  = cddy - coy                                          
 	             rsz  = cddz - coz                                          
 	             rad2 = SQRT( rsx**2 + rsy**2 + rsz**2 )

 	             if ( kk .eq. 1 ) then   ! point of reference of molec is atom closest to center.
 	              brad2 = rad2
 	             else
 	              if ( rad2 .lt. brad2 ) then
 	               brad2 = rad2
 	              endif
 	             endif
  	                    
 	    write(vvv,dp1c1) ppk,molty(ppk),kk,cox,coy,coz,rad2,brad2
 	            enddo

!-------VI.E.9. RADIAL DISTRIBUTION FUNCTION - DIST BET CLUSTER CENTER AND INDIV MOLEC
 
 	            rrb = 0d0
 	            larg= 0d0                    
 	            do rrf = 1,lrb           ! records dist of molec from center.                              
 	             rrb = rrb + bnr                                     
 	             if ( brad2 .le. rrb ) then   
 	              rgen4( rrf ) = rgen4( rrf ) + 1
 	              if    (  molty( ppk ) .eq. 1 ) then                     
 	               rwat4( rrf ) = rwat4( rrf ) + 1
 	               larg=rrb            
 	               EXIT                  !once recorded, proceeds to next molec                             
 	              elseif ( molty( ppk ) .eq. 2 ) then                   
 	               rnon4( rrf ) = rnon4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 3 ) then                     
 	               rbut4( rrf ) = rbut4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 4 ) then                   
 	               ramm4( rrf ) = ramm4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 5 ) then                   
 	               rmet4( rrf ) = rmet4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 6 ) then                   
 	               ract4( rrf ) = ract4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              elseif ( molty( ppk ) .eq. 8 ) then                   
 	               roct4( rrf ) = roct4( rrf ) + 1           
 	               larg=rrb            
 	               EXIT                                               
 	              endif                                               
 	             endif                   ! only molec within radius bin may proceed.                             
 	            enddo                    ! end of iterating radius bins.
 	 
 	            if ( larg .gt. lres ) then
 	                lres = larg
 	            endif
 	 
 	           endif                     ! only molec of current frag/cluster may proceed.
 	          enddo                      ! end of iterating all molec in current cluster.                                                     
 
!-------VI.E.10. RADIAL DISTRIBUTION FUNCTION - SUMMING COUNT FOR CLUSTERS OF SAME SIZE
 
 	          do rrf = 1,lrb                                             
 	           rmmg4( rrf ) = rmmg4( rrf ) + DBLE(rgen4( rrf ))      
 	           rmmw4( rrf ) = rmmw4( rrf ) + DBLE(rwat4( rrf ))      
 	           rmmn4( rrf ) = rmmn4( rrf ) + DBLE(rnon4( rrf ))      
 	           rmmb4( rrf ) = rmmb4( rrf ) + DBLE(rbut4( rrf ))      
 	           rmma4( rrf ) = rmma4( rrf ) + DBLE(ramm4( rrf ))      
 	           rmmm4( rrf ) = rmmm4( rrf ) + DBLE(rmet4( rrf ))      
 	           rmmc4( rrf ) = rmmc4( rrf ) + DBLE(ract4( rrf ))      
 	           rmmo4( rrf ) = rmmo4( rrf ) + DBLE(roct4( rrf ))      
 	          enddo                                                     
 	         endif                       ! only molec contained in current frag/cluster may proceed.
 	        enddo                        ! end of searching frags with current cluster size.
 	       endif                         ! only molec contained in current mer size may proceed.
 	      enddo                          ! end of iterating through all molecs.
 
!-------VI.E.11. RADIAL DISTRIBUTION FUNCTION - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
 
 	      do rrf = 1,lrb
 	       rmmg4( rrf ) = rmmg4( rrf ) / DBLE( tfrag )
 	       rmmw4( rrf ) = rmmw4( rrf ) / DBLE( tfrag )
 	       rmmn4( rrf ) = rmmn4( rrf ) / DBLE( tfrag )
 	       rmmb4( rrf ) = rmmb4( rrf ) / DBLE( tfrag )
 	       rmma4( rrf ) = rmma4( rrf ) / DBLE( tfrag )
 	       rmmm4( rrf ) = rmmm4( rrf ) / DBLE( tfrag )
 	       rmmc4( rrf ) = rmmc4( rrf ) / DBLE( tfrag )
 	       rmmo4( rrf ) = rmmo4( rrf ) / DBLE( tfrag )
 	      enddo
 
!-------VI.E.12. RDAIAL DISTRIBUTION FUNCTION - SUMMING COUNT FOR CLUSTER SIZES IN CHOSEN RANGE
 
 	      crm4( nth ) = crm4( nth ) + 1
 	      do rrf = 1,lrb
 	       genmol4( nth,rrf ) = genmol4( nth,rrf ) + rmmg4( rrf ) 
    	       sumwat4( nth,rrf ) = sumwat4( nth,rrf ) + rmmw4( rrf )     
    	       sumnon4( nth,rrf ) = sumnon4( nth,rrf ) + rmmn4( rrf )
    	       sumbut4( nth,rrf ) = sumbut4( nth,rrf ) + rmmb4( rrf )
    	       sumamm4( nth,rrf ) = sumamm4( nth,rrf ) + rmma4( rrf )
    	       summet4( nth,rrf ) = summet4( nth,rrf ) + rmmm4( rrf )
    	       sumact4( nth,rrf ) = sumact4( nth,rrf ) + rmmc4( rrf )
    	       sumoct4( nth,rrf ) = sumoct4( nth,rrf ) + rmmo4( rrf )
    	      enddo
  	     endif
    	    enddo                            ! end of iterating all cluster sizes.
  	           
 	    write( rmc,dp1h3 ) "#","BIN ID","RADIUS (A)","VOLUME (A^3)",
     &                       "N [ALL]", "g(r)",
     &                       "N [WAT]", "g(r)",
     &                       "N [NON]", "g(r)",
     &                       "N [BUT]", "g(r)",
     &                       "N [AMM]", "g(r)",
     &                       "N [MET]", "g(r)",
     &                       "N [ACT]", "g(r)",
     &                       "N [OCT]", "g(r)"

c 	    write( rmc,dfh4 ) "#"
c 	    write( rmc,dfh4 ) "#","RADIUS (nm)","DENSITY (nm^-3) =",
c     &                        "WATER","NONANE","BUTANOL","AMMONIA"
 
!-------VI.E.13. RADIAL DISTRIBUTION FUNCTION - REMOVING EXTRA VALUES	

 	    cclm3(nth) = 0
 	    do rrf = lrb,1,-1
 	     if ( genmol3(nth,rrf) .eq. 0 ) then
 	      cclm3(nth) = cclm3(nth) + 1
 	     else
 	      kin2(nth) = lrb - cclm3(nth)
 	      EXIT
 	     endif
 	    enddo

!-------VI.E.14. RADIAL DISTRIBUTION FUNCTION - AVERAGING COUNT TO REPRESENT CHOSEN RANGE
 
 	    rrb=0d0
 	    do rrf = 1,kin2(nth)
 	     rrb = rrb + bnr  
 	     if ( crm4( nth ) .eq. 0 ) then
 	      dvgen( nth,rrf ) = 0d0
 	      dvwat( nth,rrf ) = 0d0
 	      dvnon( nth,rrf ) = 0d0
 	      dvbut( nth,rrf ) = 0d0
 	      dvamm( nth,rrf ) = 0d0
 	      dvmet( nth,rrf ) = 0d0
 	      dvact( nth,rrf ) = 0d0
 	      dvoct( nth,rrf ) = 0d0
 	     else
 	      dvgen( nth,rrf ) = genmol4(nth,rrf)/DBLE(crm4(nth)) 
 	      dvwat( nth,rrf ) = sumwat4(nth,rrf)/DBLE(crm4(nth))
 	      dvnon( nth,rrf ) = sumnon4(nth,rrf)/DBLE(crm4(nth))
 	      dvbut( nth,rrf ) = sumbut4(nth,rrf)/DBLE(crm4(nth))
 	      dvamm( nth,rrf ) = sumamm4(nth,rrf)/DBLE(crm4(nth))
 	      dvmet( nth,rrf ) = summet4(nth,rrf)/DBLE(crm4(nth))
 	      dvact( nth,rrf ) = sumact4(nth,rrf)/DBLE(crm4(nth))
 	      dvoct( nth,rrf ) = sumoct4(nth,rrf)/DBLE(crm4(nth))
 	     endif
 
!-------VI.E.15. RADIAL DISTRIBUTION FUNCTION - CALCULATING VOLUME AND DENSITY
 
 	    divvol = (4d0/3d0)*3.1415927d0*(rrb*10)**3        
     &      - (4d0/3d0)*3.1415927d0*( rrb*10 - bnr*10 )**3 ! volume of bin r.
  	             
 	     ovdens = nap/((4d0/3d0)*3.1415927d0*lres**3)  ! overall density.
 	     
 	     ndgen( nth,rrf ) = dvgen( nth,rrf ) / (divvol*ovdens)
 	     ndwat( nth,rrf ) = dvwat( nth,rrf ) / (divvol*ovdens) ! distribution function.
 	     ndnon( nth,rrf ) = dvnon( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     ndbut( nth,rrf ) = dvbut( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     ndamm( nth,rrf ) = dvamm( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     ndmet( nth,rrf ) = dvmet( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     ndact( nth,rrf ) = dvact( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     ndoct( nth,rrf ) = dvoct( nth,rrf ) / (divvol*ovdens) ! distribution function. 
 	     
 	     rmd = 3066
 	     mxv = np*3
 	     rmc = rmd + mxv + nth
 	     der = rrb - bnr
 	   
 	    write(rmc,dp1c2) rrf,der*10d0,divvol,
     &                genmol4( nth,rrf ),ndgen( nth,rrf ),
     &                sumwat4( nth,rrf ),ndwat( nth,rrf ),
     &                sumnon4( nth,rrf ),ndnon( nth,rrf ),
     &                sumbut4( nth,rrf ),ndbut( nth,rrf ),
     &                sumamm4( nth,rrf ),ndamm( nth,rrf ), 
     &                summet4( nth,rrf ),ndmet( nth,rrf ), 
     &                sumact4( nth,rrf ),ndact( nth,rrf ), 
     &                sumoct4( nth,rrf ),ndoct( nth,rrf ) 

 	     rmc = rmd + mxv + nth - 2000
 	    
 	     write(rmc,*) rrf,der*10d0,divvol,ovdens,nap,lres,
     &                 genmol3( nth,rrf ),ndgen( nth,rrf ),
     &                 sumwat3( nth,rrf ),ndwat( nth,rrf ),
     &                 sumnon3( nth,rrf ),ndnon( nth,rrf ),
     &                 sumbut3( nth,rrf ),ndbut( nth,rrf ),
     &                 sumamm3( nth,rrf ),ndamm( nth,rrf ), 
     &                 summet3( nth,rrf ),ndmet( nth,rrf ), 
     &                 sumact3( nth,rrf ),ndact( nth,rrf ), 
     &                 sumoct3( nth,rrf ),ndoct( nth,rrf ) 
   	    enddo
   	   enddo
  	  endif                              ! only time matching with target time may proceed.
   	 enddo                               ! roll call of target time.
   	 endif

!-------VI.F.1. CLUSTER SIZE FREQUENCY - SETTING UP VARIABLES

    	 if ( RPSF .eq. 1 ) then    
    	 mx       = 0                        ! for growthmap output prep.
   	 my       = 0                        ! for growthmap output prep.  
   	 allfrag1 = 0                        ! collects value for cluster sizes in chosen range.
   	 plcholda = 11000                    ! placeholder.
   	 spha     = '(a1,a17,a15,a17)' 
   	 sphaa    = '(a1,a17,a22,a17)' 
 	 spca     = '(i1,f22.17,i8,i16)'     ! format specifier.
 	 spcaa    = '(i1,f22.17,f12.1,i16)'   

 	 if ( fr .eq. ffmmi ) then
 	  write( 3100,spha  ) "#","TIME (ns)","SIZE","FREQUENCY"
 	  write( 3110,spha  ) "#","TIME (ns)","SIZE","FREQUENCY"
 	  write( 3120,sphaa ) "#","TIME (ns)","SIZE (10^2)","FREQUENCY"
 	 endif

!-------VI.F.2. CLUSTER SIZE FREQUENCY - ITERATE ALL MER SIZE

  	 do mmi = 1,pop    
	  if ( tally( mmi ) .ne. 0 ) then
  	   tfrag    = tally( mmi ) / mmi  
 	   allfrag1 = allfrag1 + tfrag
 	   write( 3100,spca ) 1 , ctme , mmi , tfrag    

!-------VI.F.3. CLUSTER SIZE FREQUENCY - APPLYING PLACE HOLDERS

 	  else 
 	   tfrag  = plcholda      
 	   my     = my + 1
 	  endif
  	  write( 3110,spca ) 1 , ctme , mmi , tfrag        

!-------VI.F.4. CLUSTER SIZE FREQUENCY - PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP

 	  if ( ( mmi .eq. 1 ) .and. ( tally( mmi ) .ne. 0 ) ) then
 	   my = my + 1 
 	  endif
 	  mx  = mx + 1
 	  if ( trcnv .eq. 500 ) then 
 	   jj = 250                          ! determines cube size in growth map.
 	   pp = 3
 	  else
 	   jj = 50                           ! determines cube size in growth map.
 	   pp = 7
 	  endif 
 	  uu = (pop/jj)*2
 	  yy = DBLE(trcnv)/1000d0
 	  bq = DBLE(mx-jj)/100d0 
 	  br = DBLE(mx)/100d0
 	
 	  if ( MOD(mx, jj) .eq. 0 ) then
 	   if ( my .eq. jj ) then  
 	    write( 3120,spcaa ) 1 , ctme , bq , plcholda  
 	    write( 3120,spcaa ) 1 , ctme , br, plcholda    
 	   else       
 	    vall3 = allfrag1
  	    write( 3120,spcaa ) 1 , ctme , bq , vall3  
  	    write( 3120,spcaa ) 1 , ctme , br , vall3   
 	   endif                                                    
 	   my       = 0
 	   allfrag1 = 0
 	  endif                              ! only present mer sizes may proceed.
 	 enddo                               ! end of iterating all cluster sizes. 
 	 endif

!-------VI.G.1. MOLE FRACTION - SETTING UP VARIABLES

    	 if ( RPMOF .eq. 1 ) then    
    	 mx       = 0                        ! for growthmap output prep.             
   	 my       = 0                        ! for growthmap output prep.                   
  	 getavml1 = 0d0                      ! collects value for cluster sizes in chosen range. 
 	 getavml2 = 0d0           
 	 getavml3 = 0d0           
 	 getavml4 = 0d0           
   	 allfrag  = 0
 	 plcholdb = 2d0                      ! place holder.
 	 sphb     = '(a1,a17,a15,a23)'
 	 sphbb    = '(a1,a17,a22,a23)'
 	 spcb     = '(i1,f22.17,i8,f28.17)'  ! format specifier.      
 	 spcbb    = '(i1,f22.17,f12.1,f28.17)'      
    	        
    	 do ffi = 1,frag
    	  FRGPOOL( ffi ) = .TRUE.
    	 enddo

 	 if ( fr .eq. ffmmi ) then
 	  write( 3200,sphb ) "#","TIME (ns)","SIZE","MFRAC XWATER"
 	  write( 3201,sphb ) "#","TIME (ns)","SIZE","MFRAC XNONANE"
 	  write( 3202,sphb ) "#","TIME (ns)","SIZE","MFRAC XBUTANOL"
 	  write( 3203,sphb ) "#","TIME (ns)","SIZE","MFRAC XAMMONIA"
 	  write( 3204,sphb ) "#","TIME (ns)","SIZE","MFRAC XMETHANOL"
 	  write( 3205,sphb ) "#","TIME (ns)","SIZE","MFRAC XACETIC"
 	  write( 3206,sphb ) "#","TIME (ns)","SIZE","MFRAC XOCTANOL"
 	
 	  write( 3210,sphb ) "#","TIME (ns)","SIZE","MFRAC XWATER"
 	  write( 3211,sphb ) "#","TIME (ns)","SIZE","MFRAC XNONANE"
 	  write( 3212,sphb ) "#","TIME (ns)","SIZE","MFRAC XBUTANOL"
 	  write( 3213,sphb ) "#","TIME (ns)","SIZE","MFRAC XAMMONIA"
 	  write( 3214,sphb ) "#","TIME (ns)","SIZE","MFRAC XMETHANOL"
 	  write( 3215,sphb ) "#","TIME (ns)","SIZE","MFRAC XACETIC"
 	  write( 3216,sphb ) "#","TIME (ns)","SIZE","MFRAC XOCTANOL"
 	
 	  write( 3220,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XWATER"
 	  write( 3221,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XNONANE"
 	  write( 3222,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XBUTANOL"
 	  write( 3223,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XAMMONIA"
 	  write( 3224,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XMETHANOL"
 	  write( 3225,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XACETIC"
 	  write( 3226,sphbb)"#","TIME (ns)","SIZE(10^2)","MFRAC XOCTANOL"
 	 endif

!-------VI.G.2. MOLE FRACTION - ITERATE ALL MER SIZE

  	 do mmi = 1,pop    
	  if ( tally( mmi ) .ne. 0 ) then
   	   tmlfrcn1 = 0d0                    ! collects value for clusters of same size.   
   	   tmlfrcn2 = 0d0          
   	   tmlfrcn3 = 0d0          
   	   tmlfrcn4 = 0d0          
   	   tmlfrcn5 = 0d0          
   	   tmlfrcn6 = 0d0          
   	   tmlfrcn8 = 0d0          
  	   tfrag = tally( mmi )/mmi  

!-------VI.G.3. MOLE FRACTION - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE

   	   do ppi = 1,pop  
	    if ( mer( ppi ) .eq. mmi ) then
	     do ffj = 1,frag  
	      if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &             ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
	
	       FRGPOOL( ffj ) = .FALSE.  
	       typcnt1 = 0              
	       typcnt2 = 0             
	       typcnt3 = 0             
	       typcnt4 = 0               
	       typcnt5 = 0               
	       typcnt6 = 0               
	       typcnt8 = 0               
                         
!-------VI.G.4. MOLE FRACTION - INDEXING MOLECULES INVOLVED IN LOCATED FRAG

 	       do ppj = 1,pop                                             
 	        if ( nuc( ppj ) .eq. ffj ) then
 	         if     ( molty( ppj ) .eq. 1 ) then
 	          typcnt1 = typcnt1 + 1
 	         elseif ( molty( ppj ) .eq. 2 ) then
 	          typcnt2 = typcnt2 + 1
 	         elseif ( molty( ppj ) .eq. 3 ) then
 	          typcnt3 = typcnt3 + 1
 	         elseif ( molty( ppj ) .eq. 4 ) then
 	          typcnt4 = typcnt4 + 1
 	         elseif ( molty( ppj ) .eq. 5 ) then
 	          typcnt5 = typcnt5 + 1
 	         elseif ( molty( ppj ) .eq. 6 ) then
 	          typcnt6 = typcnt6 + 1
 	         elseif ( molty( ppj ) .eq. 8 ) then
 	          typcnt8 = typcnt8 + 1
 	         endif
 	        endif                                
 	       enddo                                   

!-------VI.G.5. MOLE FRACTION - CALCULATING MOLE FRACTION OF CURRENT FRAG

   	       mltot=typcnt1+typcnt2+typcnt3+typcnt4+typcnt5+typcnt6
   	
   	       mlfrcn1 = DBLE( typcnt1 ) / DBLE( mltot )
   	       mlfrcn2 = DBLE( typcnt2 ) / DBLE( mltot )
   	       mlfrcn3 = DBLE( typcnt3 ) / DBLE( mltot )
   	       mlfrcn4 = DBLE( typcnt4 ) / DBLE( mltot )
   	       mlfrcn5 = DBLE( typcnt5 ) / DBLE( mltot )
   	       mlfrcn6 = DBLE( typcnt6 ) / DBLE( mltot )
   	       mlfrcn8 = DBLE( typcnt8 ) / DBLE( mltot )

!-------VI.G.6. MOLE FRACTION - SUMMING COUNT FOR CLUSTERS OF SAME SIZE

   	       tmlfrcn1 = tmlfrcn1 + mlfrcn1    
   	       tmlfrcn2 = tmlfrcn2 + mlfrcn2    
   	       tmlfrcn3 = tmlfrcn3 + mlfrcn3    
   	       tmlfrcn4 = tmlfrcn4 + mlfrcn4    
   	       tmlfrcn5 = tmlfrcn5 + mlfrcn5    
   	       tmlfrcn6 = tmlfrcn6 + mlfrcn6    
   	       tmlfrcn8 = tmlfrcn8 + mlfrcn8    
   	      endif                          ! only molec in current frag/cluster may proceed.
   	      enddo                          ! end of searching frag with current cluster size.
   	    endif                            ! only molec contained in current mer size may proceed.
 	   enddo                             ! end of iterating all molec in current cluster.

!-------VI.G.7. MOLE FRACTION - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE

 	   avrgml1 = tmlfrcn1 / DBLE( tfrag )    
 	   avrgml2 = tmlfrcn2 / DBLE( tfrag )    
 	   avrgml3 = tmlfrcn3 / DBLE( tfrag )    
 	   avrgml4 = tmlfrcn4 / DBLE( tfrag )    
 	   avrgml5 = tmlfrcn5 / DBLE( tfrag )    
 	   avrgml6 = tmlfrcn6 / DBLE( tfrag )    
 	   avrgml8 = tmlfrcn8 / DBLE( tfrag )    
 	
 	   write( 3200, spcb ) 1 , ctme , mmi , avrgml1  
 	   write( 3201, spcb ) 1 , ctme , mmi , avrgml2  
 	   write( 3202, spcb ) 1 , ctme , mmi , avrgml3  
 	   write( 3203, spcb ) 1 , ctme , mmi , avrgml4  
 	   write( 3204, spcb ) 1 , ctme , mmi , avrgml5  
 	   write( 3205, spcb ) 1 , ctme , mmi , avrgml6  
 	   write( 3206, spcb ) 1 , ctme , mmi , avrgml8  

!-------VI.G.8. MOLE FRACTION - SUMMING COUNT TO REPRESENT CHOSEN RANGE

 	   if ( mmi .ne. 1 ) then
 	    getavml1 = getavml1 + avrgml1 * DBLE( tfrag )       
 	    getavml2 = getavml2 + avrgml2 * DBLE( tfrag )       
 	    getavml3 = getavml3 + avrgml3 * DBLE( tfrag )       
 	    getavml4 = getavml4 + avrgml4 * DBLE( tfrag )       
 	    getavml5 = getavml5 + avrgml5 * DBLE( tfrag )       
 	    getavml6 = getavml6 + avrgml6 * DBLE( tfrag )       
 	    getavml8 = getavml8 + avrgml8 * DBLE( tfrag )       
 	    allfrag  = allfrag  + tfrag                       
 	   endif

!-------VI.G.9. MOLE FRACTION - APPLYING PLACE HOLDERS

 	  else 
 	   avrgml1 = plcholdb        
 	   avrgml2 = plcholdb        
 	   avrgml3 = plcholdb        
 	   avrgml4 = plcholdb        
 	   avrgml5 = plcholdb        
 	   avrgml6 = plcholdb        
 	   avrgml8 = plcholdb        
 	   my      = my + 1
 	  endif
 	  
  	  write( 3210,spcb ) 1 , ctme , mmi , avrgml1      
  	  write( 3211,spcb ) 1 , ctme , mmi , avrgml2      
  	  write( 3212,spcb ) 1 , ctme , mmi , avrgml3      
  	  write( 3213,spcb ) 1 , ctme , mmi , avrgml4      
  	  write( 3214,spcb ) 1 , ctme , mmi , avrgml5     
  	  write( 3215,spcb ) 1 , ctme , mmi , avrgml6     
  	  write( 3216,spcb ) 1 , ctme , mmi , avrgml8     

!-------VI.G.10. MOLE FRACTION - PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP

 	  if ( ( mmi .eq. 1 ) .and. ( tally( mmi ) .ne. 0 ) ) then
 	   my = my + 1 
 	  endif
 	  mx = mx + 1
 	  if ( trcnv .eq. 500 ) then 
 	   jj = 250                          ! determines cube size in growth map.
 	   pp = 3                               
 	  else                                  
 	   jj = 50                           ! determines cube size in growth map.
 	   pp = 7
 	  endif 
 	  uu = (pop/jj)*2
 	  yy = DBLE(trcnv)/1000d0
 	  bq = DBLE(mx-jj)/100d0 
 	  br = DBLE(mx)/100d0
 	
 	  if ( MOD(mx, jj) .eq. 0 ) then
 	   if ( my .eq. jj ) then  
 	    write( 3220,spcbb ) 1 , ctme , bq , 2d0      
 	    write( 3220,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3221,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3221,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3222,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3222,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3223,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3223,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3224,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3224,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3225,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3225,spcbb ) 1 , ctme , br , 2d0     
 	    write( 3226,spcbb ) 1 , ctme , bq , 2d0     
 	    write( 3226,spcbb ) 1 , ctme , br , 2d0     
 	   else       
 	    vall2W = getavml1 / DBLE( allfrag )
 	    vall2N = getavml2 / DBLE( allfrag )
 	    vall2B = getavml3 / DBLE( allfrag )
 	    vall2A = getavml4 / DBLE( allfrag )
 	    vall2M = getavml5 / DBLE( allfrag )
 	    vall2C = getavml6 / DBLE( allfrag )
 	    vall2O = getavml8 / DBLE( allfrag )

  	    write( 3220,spcbb ) 1 , ctme , bq , vall2W  
  	    write( 3220,spcbb ) 1 , ctme , br , vall2W  
  	    write( 3221,spcbb ) 1 , ctme , bq , vall2N  
  	    write( 3221,spcbb ) 1 , ctme , br , vall2N  
  	    write( 3222,spcbb ) 1 , ctme , bq , vall2B  
  	    write( 3222,spcbb ) 1 , ctme , br , vall2B  
  	    write( 3223,spcbb ) 1 , ctme , bq , vall2A  
  	    write( 3223,spcbb ) 1 , ctme , br , vall2A  
  	    write( 3224,spcbb ) 1 , ctme , bq , vall2M  
  	    write( 3224,spcbb ) 1 , ctme , br , vall2M  
  	    write( 3225,spcbb ) 1 , ctme , bq , vall2C  
  	    write( 3225,spcbb ) 1 , ctme , br , vall2C  
  	    write( 3226,spcbb ) 1 , ctme , bq , vall2O  
  	    write( 3226,spcbb ) 1 , ctme , br , vall2O  
 	   endif                                                    

 	   my = 0
 	   getavml1 = 0d0   
 	   getavml2 = 0d0   
 	   getavml3 = 0d0   
 	   getavml4 = 0d0   
 	   getavml5 = 0d0   
 	   getavml6 = 0d0   
 	   getavml8 = 0d0   
 	   allfrag  = 0     
 	  endif                      
    	 enddo                               ! end of iterating all cluster sizes.
    	 endif

!-------VI.H.1. CLUSTER MASS - SETTING UP VARIABLES

	 if ( RPM .eq. 1  ) then    
   	 l        = 0
   	 mq       = 0                   
    	 mx       = 0                        ! for growthmap output prep.  
   	 my       = 0                        ! for growthmap output prep.        
   	 getavmas = 0d0                      ! collects value for cluster sizes in chosen range. 
   	 getavupm = 0d0          
  	 getavlom = 0d0     
   	 allfrag  = 0             
  	 plcholdc = 500000d0                    ! placeholder.
 	 sphc     = '(a1,a17,a15,a40)'
 	 sphcc    = '(a1,a17,a22,a40)'
  	 spcc     = '(i1,f22.17,i8,f50.36)'     ! format specifier.
  	 spccc    = '(i1,f22.17,f12.1,f50.36)'  ! format specifier.
  	 avog     = 6.0221409d0*10d0**23        ! avogadro's number.
    	
    	 do ffi = 1,frag
    	  FRGPOOL( ffi ) = .TRUE.
    	 enddo

   	 do ppi = 1,pop                      ! determines smaller and larger molar mass.
    	  mssm = 0d0
   	  if     ( ( molty( ppi ) .eq. 1 )  .and. 
     &             ( l            .ne. 4 ) ) then
   	   l  = 4
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 2 )  .and. 
     &             ( l            .ne. 9 ) ) then
   	   l  = 9
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 3 )  .and. 
     &             ( l            .ne. 6 ) ) then
   	   l  = 6
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 4 )  .and. 
     &             ( l            .ne. 5 ) ) then
   	   l  = 5
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 5 )  .and. 
     &             ( l            .ne. 3 ) ) then
   	   l  = 3
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 6 )  .and. 
     &             ( l            .ne. 5 ) ) then
   	   l  = 5
   	   mq = mq + 1
   	  elseif ( ( molty( ppi ) .eq. 8 )  .and. 
     &             ( l            .ne. 10 ) ) then
   	   l  = 10
   	   mq = mq + 1
   	  endif
   	
   	  do k = 1,l                                              
   	   mssm = mssm + mss( molty( ppi ),k )  
   	  enddo
   	  
   	  if ( mq .eq. 1 ) then
   	   massl = mssm
   	   massh = mssm
   	  else 
   	   if ( mssm .lt. massl ) then
   	    massl = mssm                     ! determine smaller molar mass
   	   else
   	    massh = mssm                     ! determine larger molar mass.
   	   endif
   	   EXIT                              ! done determining contained species.
   	  endif
   	 enddo
         
 	 if ( fr .eq. ffmmi ) then
 	  write( 3300,sphc ) "#","TIME (ns)","SIZE","MASS (g)"
 	  write( 3301,sphc ) "#","TIME (ns)","SIZE",
     &                       "MASS FRACTION (g/max g)"
 	  write( 3302,sphc ) "#","TIME (ns)","SIZE",
     &                       "MASS FRACTION (g/min g)"
 	
 	  write( 3310,sphc ) "#","TIME (ns)","SIZE","MASS (g)"
 	  write( 3311,sphc ) "#","TIME (ns)","SIZE",
     &                       "MASS FRACTION (g/max g)"
 	  write( 3312,sphc ) "#","TIME (ns)","SIZE",
     &                       "MASS FRACTION (g/min g)"
 	
 	  write( 3320,sphcc ) "#","TIME (ns)","SIZE (10^2)",
     &                        "MASS (10^-20 g)"
 	  write( 3321,sphcc ) "#","TIME (ns)","SIZE (10^2)",
     &                        "MASS FRACTION (g/max g)"
 	  write( 3322,sphcc ) "#","TIME (ns)","SIZE (10^2)",
     &                        "MASS FRACTION (g/min g)"
 	 endif

!-------VI.H.2. CLUSTER MASS - ITERATE ALL MER SIZE

  	 do mmi = 1,pop    
	  if ( tally( mmi ) .ne. 0 ) then
  	   tfrag = tally( mmi )/mmi  
   	   tmas = 0d0                        ! collects value for clusters of same size.
  	   tupl = 0d0              
	   tlol = 0d0              
   	   uprlim = massh*DBLE( mmi )
  	   lorlim = massl*DBLE( mmi )        

!-------VI.H.3 CLUSTER MASS - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE

   	   do ppi = 1,pop  
	    if ( mer( ppi ) .eq. mmi ) then
	     do ffj = 1,frag  
	      if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &             ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
	
	       FRGPOOL( ffj ) = .FALSE.  
	       mssm = 0d0

!-------VI.H.4 CLUSTER MASS - INDEXING MOLECULES INVOLVED IN LOCATED FRAG

   	       do ppj = 1,pop                                        
   	        if ( nuc( ppj ) .eq. ffj ) then
   	         if ( molty( ppj ) .eq. 1 ) then
   	          l = 4
   	         elseif ( molty( ppj ) .eq. 2 ) then
   	          l = 9
   	         elseif ( molty( ppj ) .eq. 3 ) then
   	          l = 6
   	         elseif ( molty( ppj ) .eq. 4 ) then
   	          l = 5
   	         elseif ( molty( ppj ) .eq. 5 ) then
   	          l = 3
   	         elseif ( molty( ppj ) .eq. 6 ) then
   	          l = 5
   	         elseif ( molty( ppj ) .eq. 8 ) then
   	          l = 10
   	         endif

!-------VI.H.5 CLUSTER MASS - CALCULATING MASS OF CURRENT FRAG 

   	         do k = 1,l                                              
   	          mssm = mssm + mss( molty( ppj ),k )                       
   	         enddo                                                
   	        endif
   	       enddo
   	       upl = mssm / DBLE( uprlim )   ! determines mass fraction (over upper limit).
   	       lol = mssm / DBLE( lorlim )   ! determines mass fraction (over lower limit).  

!-------VI.H.6. CLUSTER MASS - SUMMING COUNT FOR CLUSTERS OF SAME SIZE

   	       tmas = tmas + mssm 
   	       tupl = tupl + upl
   	       tlol = tlol + lol
   	      endif                          ! only molec in current frag/cluster may proceed.
   	     enddo                           ! end of searching frag with current cluster size.
   	    endif                            ! only molec contained in current mer size may proceed.
   	   enddo                             ! end of iterating all molec in current cluster.

!-------VI.H.7. CLUSTER MASS - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE

 	   avrgms = tmas   / DBLE( tfrag )    
 	   avrgup = tupl   / DBLE( tfrag )    
 	   avrglo = tlol   / DBLE( tfrag )    
 	   relmas = avrgms / avog
 	
 	   write( 3300,spcc ) 1 , ctme , mmi , relmas   
 	   write( 3301,spcc ) 1 , ctme , mmi , avrgup   
 	   write( 3302,spcc ) 1 , ctme , mmi , avrglo   

!-------VI.H.8. CLUSTER MASS - SUMMING COUNT TO REPRESENT CHOSEN RANGE

 	   if ( mmi .ne. 1 ) then
 	    getavmas = getavmas + avrgms * DBLE( tfrag )       
 	    getavupm = getavupm + avrgup * DBLE( tfrag )       
 	    getavlom = getavlom + avrglo * DBLE( tfrag )       
 	    allfrag  = allfrag  + tfrag                       
 	   endif

!-------VI.H.9. CLUSTER MASS - APPLYING PLACE HOLDERS

 	  else 
 	   avrgms = plcholdc 
 	   avrgup = plcholdc
 	   avrglo = plcholdc
 	   relmas = avrgms  
 	   my     = my + 1
 	  endif
 	
  	  write( 3310,spcc ) 1 , ctme , mmi , relmas     
  	  write( 3311,spcc ) 1 , ctme , mmi , avrgup      
  	  write( 3312,spcc ) 1 , ctme , mmi , avrglo      

!-------VI.H.10. CLUSTER MASS - PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP

 	  if ( ( mmi .eq. 1 ) .and. ( tally( mmi ) .ne. 0 ) ) then
 	   my = my + 1 
 	  endif
 	  mx = mx + 1
 	  if ( trcnv .eq. 500 ) then 
 	   jj = 250                          ! determines cube size in growth map.
 	   pp = 3
 	  else
 	   jj = 50                           ! determines cube size in growth map.
 	   pp = 7
 	  endif 
 	  uu = (pop/jj)*2
 	  yy = DBLE(trcnv)/1000d0
 	  bq = DBLE(mx-jj)/100d0 
 	  br = DBLE(mx)/100d0
 	
 	  if ( MOD(mx, jj) .eq. 0 ) then
 	   if ( my .eq. jj ) then  
 	    write( 3320,spccc ) 1 , ctme , bq , plcholdc 
 	    write( 3320,spccc ) 1 , ctme , br , plcholdc 
 	    write( 3321,spccc ) 1 , ctme , bq , plcholdc    
 	    write( 3321,spccc ) 1 , ctme , br , plcholdc    
 	    write( 3322,spccc ) 1 , ctme , bq , plcholdc    
 	    write( 3322,spccc ) 1 , ctme , br , plcholdc    
 	   else       
 	    vall4 = getavmas / DBLE( allfrag )
 	    vall5 = vall4    / 1000d0    !avog
 	    vall6 = getavupm / DBLE( allfrag )
 	    vall7 = getavlom / DBLE( allfrag )

  	    write( 3320,spccc ) 1 , ctme , bq , vall5  
  	    write( 3320,spccc ) 1 , ctme , br , vall5  
  	    write( 3321,spccc ) 1 , ctme , bq , vall6  
  	    write( 3321,spccc ) 1 , ctme , br , vall6  
  	    write( 3322,spccc ) 1 , ctme , bq , vall7  
  	    write( 3322,spccc ) 1 , ctme , br , vall7  
 	   endif                                                    

 	   my       = 0
 	   allfrag  = 0     
 	   getavmas = 0d0   
 	   getavupm = 0d0   
 	   getavlom = 0d0   
 	  endif
 	 enddo                               ! end of iterating all cluster sizes.
 	 endif

!-------VI.I.1. CLUSTER RADIUS - SETTING UP VARIABLES

	 if ( RPR .eq. 1  ) then    
    	 mx       = 0                        ! for growthmap output prep. 
   	 my       = 0                        ! for growthmap output prep. 
   	 getavrad = 0d0                      ! collects value for cluster sizes in chosen range. 
   	 allfrag  = 0             
    	 plcholdd = 200d0                    ! placeholder.
 	 sphd     = '(a1,a17,a15,a25)'
 	 sphdd    = '(a1,a17,a22,a25)'
    	 spcd     = '(i1,f22.17,i8,f29.17)'  ! format specifier.
    	 spcdd    = '(i1,f22.17,f12.1,f29.17)' 
    	 
    	 do ffi = 1,frag
    	  FRGPOOL( ffi ) = .TRUE.
    	 enddo

 	 if ( fr .eq. ffmmi ) then
 	  write( 3400,sphd  ) "#","TIME (ns)","SIZE","RADIUS (nm)"
 	  write( 3410,sphd  ) "#","TIME (ns)","SIZE","RADIUS (nm)"
 	  write( 3420,sphdd ) "#","TIME (ns)","SIZE (10^2)","RADIUS (nm)"
 	 endif
 
!-------VI.I.2. CLUSTER RADIUS - ITERATE ALL MER SIZE

 	 do mmi = 1,pop    
 	  if ( tally( mmi ) .ne. 0 ) then
 	   tfrag = tally( mmi )/mmi  
 	   trad = 0d0                        ! collects value for clusters of same size.
 
!-------VI.I.3. CLUSTER RADIUS - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE

   	   do ppi = 1,pop  
 	    if ( mer( ppi ) .eq. mmi ) then
 	     do ffj = 1,frag  
 	      if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &             ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	
 	       FRGPOOL( ffj ) = .FALSE.  
 	       mssx = 0d0              
 	       mssy = 0d0              
 	       mssz = 0d0              
 	       mssm = 0d0
                         
!-------VI.I.4. CLUSTER RADIUS - INDEXING MOLECULES INVOLVED IN LOCATED FRAG

 	       do ppj = 1,pop                                            
 	        if ( nuc( ppj ) .eq. ffj ) then
 	         if ( molty( ppj ) .eq. 1 ) then
 	           l = 4
 	         elseif ( molty( ppj ) .eq. 2 ) then
 	           l = 9
 	         elseif ( molty( ppj ) .eq. 3 ) then
 	           l = 6
 	         elseif ( molty( ppj ) .eq. 4 ) then
 	           l = 5
 	         elseif ( molty( ppj ) .eq. 5 ) then
 	           l = 3
 	         elseif ( molty( ppj ) .eq. 6 ) then
 	           l = 5
 	         elseif ( molty( ppj ) .eq. 8 ) then
 	           l = 10
 	         endif

!-------VI.I.5. CLUSTER RADIUS - CALCULATING CENTER OF MASS OF CURRENT FRAG
 
 	         do k = 1,l                                              
 	          mssx = mssx + mss( molty(ppj),k )*x( ppj,molty(ppj),k ) 
 	          mssy = mssy + mss( molty(ppj),k )*y( ppj,molty(ppj),k ) 
 	          mssz = mssz + mss( molty(ppj),k )*z( ppj,molty(ppj),k ) 
 	          mssm = mssm + mss( molty(ppj),k )                          
 	         enddo                                                
 	        endif                                                 
 	       enddo
 	       cmmx = mssx / mssm                                  
 	       cmmy = mssy / mssm                                          
 	       cmmz = mssz / mssm                                          
 
!-------VI.I.6. CLUSTER RADIUS - CALCULATING CENTER OF MASS OF INDIV MOLEC
 
  	       tmoi = 0d0                                                 
  	       do ppk = 1,pop                                            
   	        msx = 0d0                                                
   	        msy = 0d0                                                
   	        msz = 0d0                                                
   	        msm = 0d0                                                
 	        if ( nuc( ppk ) .eq. ffj ) then
 	         if ( molty( ppk ) .eq. 1 ) then
 	          li = 4
 	         elseif ( molty(ppk) .eq. 2 ) then
 	          li = 9
 	         elseif ( molty(ppk) .eq. 3 ) then
 	          li = 6
 	         elseif ( molty(ppk) .eq. 4 ) then
 	          li = 5
 	         elseif ( molty(ppk) .eq. 5 ) then
 	          li = 3
 	         elseif ( molty(ppk) .eq. 6 ) then
 	          li = 5
 	         elseif ( molty(ppk) .eq. 8 ) then
 	          li = 10
 	         endif
 	
 	         do kk = 1,li
 	          msx = msx + mss( molty(ppk),kk )*x( ppk,molty(ppk),kk )   
 	          msy = msy + mss( molty(ppk),kk )*y( ppk,molty(ppk),kk )   
 	          msz = msz + mss( molty(ppk),kk )*z( ppk,molty(ppk),kk )   
 	          msm = msm + mss( molty(ppk),kk )                           
 	         enddo
  	         cmx = msx / msm                                           
   	         cmy = msy / msm                                           
 	         cmz = msz / msm                                           
 
!-------VI.I.7. CLUSTER RADIUS - CALCULATING RADIUS FROM MOMENT OF INERTIA
 
 	         rrx  = cmmx - cmx           ! needed to get moi of indiv molec.                            
  	         rry  = cmmy - cmy                                      
 	         rrz  = cmmz - cmz                                         
 	         rad2 = rrx**2 + rry**2 + rrz**2                             
  	         moi  = rad2 * msm                                         
  	         tmoi = tmoi + moi           ! determines moi of current frag/cluster.                         
  	        endif                                                   
  	       enddo                                                     
  	       rad = ( tmoi/mssm )**.5d0     ! determines radius of current frag/cluster.                             
               
!-------VI.I.8. CLUSTER RADIUS - SUMMING COUNT FOR CLUSTERS OF SAME SIZE

  	       trad = trad + rad                                           
   	      endif                          ! only molec in current frag/cluster may proceed.
   	     enddo                           ! end of searching frag with current cluster size.
   	    endif                            ! only molec contained in current mer size may proceed.
   	   enddo                             ! end of iterating all molec in current cluster.
 
!-------VI.I.9. CLUSTER RADIUS - AVERAGING COUNT TO REPRESENT A CLUSTER SIZE
 
 	   avrgrd  = trad / DBLE( tfrag )    
  	
 	   write( 3400,spcd ) 1 , ctme , mmi , avrgrd   

!-------VI.I.10. CLUSTER RADIUS - SUMMING COUNT TO REPRESENT CHOSEN RANGE 
 
 	   if ( mmi .ne. 1 ) then
 	    allfrag  = allfrag  + tfrag                       
 	    getavrad = getavrad + avrgrd * DBLE( tfrag )       
 	   endif

!-------VI.I.11. CLUSTER RADIUS - APPLYING PLACEHOLDERS

 	  else 
 	   avrgrd  = plcholdd
 	   my = my + 1
 	  endif
  	
  	  write( 3410,spcd ) 1 , ctme , mmi , avrgrd       

!-------VI.I.12. CLUSTER RADIUS - PREPARING COMPATIBLE OUTPUT FOR GROWTH MAP

 	  if ( ( mmi .eq. 1 ) .and. ( tally( mmi ) .ne. 0 ) ) then
 	   my = my + 1 
 	  endif
 	  mx = mx + 1
 	  if ( trcnv .eq. 500 ) then 
 	   jj = 250                            ! determine cube size in growth map. 
 	   pp = 3
 	  else
 	   jj = 50                             ! determine cube size in growth map.
 	   pp = 7
 	  endif 
 	  uu = (pop/jj)*2
 	  yy = DBLE(trcnv)/1000d0
 	  bq = DBLE(mx-jj)/100d0 
 	  br = DBLE(mx)/100d0
 	
 	  if ( MOD( mx, jj ) .eq. 0 ) then
 	   if ( my .eq. jj ) then  
 	    write( 3420,spcdd ) 1 , ctme , bq , plcholdd    
 	    write( 3420,spcdd ) 1 , ctme , br , plcholdd  
 	   else       
 	    vall1 = getavrad / DBLE( allfrag )
 
  	    write( 3420,spcdd ) 1 , ctme , bq , vall1
  	    write( 3420,spcdd ) 1 , ctme , br , vall1  
 	   endif                                                    

 	   my = 0
 	   allfrag  = 0     
 	   getavrad = 0d0   
 	  endif
 	 enddo                                 ! end of iterating all cluster sizes.
 	 endif

!-------VI.J.1. NUCLEATION RATE - SETTING UP VARIABLES 
 
  	 if ( RPNR .eq. 1 ) then
    	 sphe = '(a1,a16,a30,i3)'
   	 spce = '(f22.17,i15)'
 	 bbx = 3500
 	
 	 do mp = 1,18
 	  cth = 5*mp
 	  bbx = bbx + 5
 	  cee( cth ) = 0
 	  if ( fr .eq. ffmmi ) then
 	   write( bbx,sphe ) "#","TIME (ns)","NO. OF CLUSTERS > N =",cth
 	  endif
 	 enddo
  	
    	 do ffi = 1,frag
    	  FRGPOOL( ffi ) = .TRUE.
    	 enddo
         
!-------VI.J.2. NUCLEATION RATE - ITERATE ALL MER SIZE

  	 do mmi = 5,pop    
  	  if ( tally( mmi ) .ne. 0 ) then

!-------VI.J.3. NUCLEATION RATE - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE

   	   do ppi=1,pop  
   	    if ( mer(ppi) .eq. mmi ) then
   	     do ffj=1,frag  
   	      if ( ( nuc(ppi) .eq. ffj ) .and.                
     &             ( FRGPOOL(ffj) .eqv. .TRUE. ) ) then
   	
   	       FRGPOOL(ffj)=.FALSE.  
        
!-------VI.J.4. NUCLEATION RATE - ACCOUNTING FRAG IN THRESHOLD BIN

 	       do mp = 1,18
 	        cth = 5*mp
 	        if ( mmi .gt. cth ) then
 	
 	         cee( cth ) = cee( cth ) + 1
 	        endif
 	       enddo
 	
 	      endif                            ! only molec in current fag/cluster may proceed.
 	     enddo                             ! end of searching frag with current cluster size.
 	    endif                              ! only molec contained in current mer size may proceed.
 	   enddo                               ! end of iterating all molec in current cluster.
 	   
 	  endif
 	 enddo                                 ! end of iterating all cluster sizes.
 	
  	 bbx = 3500 
  	 do mp = 1,18
  	  cth = 5*mp
  	  bbx = bbx + 5
  	  write( bbx,spce ) ctme,cee( cth ) 
  	 enddo
  	 endif

!-------VI.K.1. CRITICAL CLUSTER SIZE - SETTING UP VARIABLES
       
   	 if ( RPCS .eq. 1 ) then 
   	 manc = 80       !maximum k change considered       
   	 tidv = 50d0     !each frame .eq. 50 ps 
   	 csh1 = '(a1,a6,a28,a30,a34)'
   	 csc1 = '(i6,f20.4,f30.4,f29.4)'
   	
   	 do ffi = 1,frag
   	  FRGPOOL( ffi ) = .TRUE.
   	 enddo

   	 if ( fr .eq. ffmmi ) then
   	  write(3600,csh1) "#","SIZE","MEAN DECAY RATE (ps^-1)",
     &                                "MEAN GROWTH RATE (ps^-1)", 
     &                                "MEAN SIZE CHANGE RATE (ps^-1)"
   	  do binc = 1,manc                
   	    decaytally( mmi,binc ) = 0
   	   growthtally( mmi,binc ) = 0
   	  enddo
   	  zerotally( mmi ) = 0

!-------VI.K.2. CRITICAL CLUSTER SIZE - ITERATE MER SIZE

   	 else
   	  do mmi = 2,80  
   	   if ( prevtally( mmi ) .ne. 0 ) then 

!-------VI.K.3. CRITICAL CLUSTER SIZE - SEARCHING FRAGS INVOLVED IN CURRENT SIZE

   	    do ppi = 1,pop  
   	     if ( prevmer( ppi ) .eq. mmi ) then
 	      do ffj = 1,prevfrag  
 	       if ( prevnuc( ppi ) .eq. ffj ) then
 	   
   	         decaysum( ppi ) = 0
   	        growthsum( ppi ) = 0

!-------VI.K.4. CRITICAL CLUSTER SIZE - SEARCHING DECAY AND GROWTH OCCURRENCE

	        do ppj = 1,pop                                               
	         if ( ( prevnuc(ppj) .eq. ffj  ) .and.       
     &                ( nuc(ppi) .ne. nuc(ppj) ) ) then
	           decaysum(ppi) = decaysum(ppi)  + 1
	         elseif ( ( prevnuc(ppj) .ne. ffj)  .and.
     &                    ( nuc(ppi) .eq. nuc(ppj) ) ) then
	          growthsum(ppi) = growthsum(ppi) + 1
	         endif
    	        enddo
   	       endif                       ! only molec contained in current frag/cluster may proceed.
   	      enddo                        ! end of searching frags with current cluster size.
   	     endif                         ! only molec contained in current mer size may proceed.
   	    enddo                          ! end of iterating through all molecs.   

!-------VI.K.5. CRITICAL CLUSTER SIZE - DETERMINE ACCEPTABLE COUNT

   	    do ffj = 1,prevfrag  
   	     fill = 0
   	     pas  = 0
   	     do ppi = 1,pop  
   	      if ( ( prevnuc( ppi ) .eq. ffj ) .and.
     &             ( prevmer( ppi ) .eq. mmi ) ) then
   	       fill = fill + 1
   	       if ( pas .eq. 0 ) then
   	        pas = 1
   	        decaynuc(ffj)  = decaysum(ppi)
   	        growthnuc(ffj) = growthsum(ppi)
   	       elseif (decaysum(ppi) .lt. decaynuc(ffj) ) then
   	        decaynuc(ffj) = decaysum(ppi)
   	        growthnuc(ffj) = growthsum(ppi)
   	       elseif ( (decaysum(ppi)  .eq. decaynuc(ffj)  ) .and.
     &                  (growthsum(ppi) .lt. growthnuc(ffj) ) ) then
   	        growthnuc(ffj)=growthsum(ppi)
   	       endif
   	       if (fill .eq. mmi) then
   	        EXIT
   	       endif
   	      endif
   	     enddo

!-------VI.K.6. CRITICAL CLUSTER SIZE - TALLY NO CHANGE FOR EACH MER SIZE

   	     do ppi = 1,pop  
   	      if ( ( prevnuc( ppi ) .eq. ffj ) .and.
     &             ( prevmer( ppi ) .eq. mmi ) .and. 
     &             ( FRGPOOL( ffj ) .eqv. .TRUE.) ) then
   	       FRGPOOL( ffj ) = .FALSE.
   	       if ( (  decaynuc(ffj) .eq. 0 ) .and. 
     &              ( growthnuc(ffj) .eq. 0 ) ) then
   	        zerotally(mmi) = zerotally(mmi) + 1

!-------VI.K.7. CRITICAL CLUSTER SIZE - TALLY CHANGE COUNT FOR EACH MER SIZE

   	       else
   	
   	        do binc = 1,manc
   	         if ( decaynuc(ffj) .eq. binc) then
   	          decaytally(mmi,binc)  = decaytally(mmi,binc)  + 1
   	          EXIT
   	         endif
   	        enddo
   	
   	        do binc = 1,manc
   	         if ( growthnuc(ffj) .eq. binc) then
   	          growthtally(mmi,binc) = growthtally(mmi,binc) + 1
   	          EXIT
   	         endif
   	        enddo
   	       endif
   	      endif
   	     enddo
   	    enddo
   	   endif
   	  enddo                            ! end of iterating all cluster sizes.
   	 endif

!-------VI.K.8. CRITICAL CLUSTER SIZE - TIME AVERAGE FOR EACH MER SIZE

   	 if ( fr .eq. 100 ) then
   	  do mmi = 2,80          
   	   bankgrowth  = 0
   	   bankdecay   = 0
   	   changesum   = 0
   	   derate(mmi) = 0
   	   grrate(mmi) = 0 
   	  
   	   do binc = 1,manc
   	    avgdty(mmi,binc) = DBLE(decaytally(mmi,binc))/100d0         ! 100 frames.
   	    avggty(mmi,binc) = DBLE(growthtally(mmi,binc))/100d0
   	    bankgrowth = bankgrowth + avggty(mmi,binc) 
   	    bankdecay  = bankdecay  + avgdty(mmi,binc)
   	   enddo
   	   avgzty(mmi) = zerotally(mmi)/100d0
   	   bankzero = avgzty(mmi)
   	   changesum = bankdecay + bankgrowth + bankzero

!-------VI.K.9. CRITICAL CLUSTER SIZE - TRANSITION PROBABILITY AND CHANGE RATE

   	   do binc = 1,manc
   	    tdprob(mmi,binc) = ( avgdty(mmi,binc) )/changesum
   	    tgprob(mmi,binc) = ( avggty(mmi,binc) )/changesum
   	    derate(mmi) = derate(mmi) + DBLE(binc)*tdprob(mmi,binc)
   	    grrate(mmi) = grrate(mmi) + DBLE(binc)*tgprob(mmi,binc)
   	   enddo
   	    
   	   derate(mmi) = derate(mmi)/tidv    ! unit becomes (ps^-1), 
   	   grrate(mmi) = grrate(mmi)/tidv    ! instead of (50ps^-1).
   	   sumchangerate(mmi) = -derate(mmi) + grrate(mmi)   
   	   write(3600,csc1) mmi,-derate(mmi),                  
     &                           grrate(mmi),sumchangerate(mmi)
   	  enddo
   	 endif

!-------VI.K.10. CRITICAL CLUSTER SIZE - RECORD MER, FRAG, AND TALLY ID FOR NEXT TIME STEP

   	 do ppi = 1,pop
   	  prevnuc(ppi)=nuc(ppi)
   	  prevmer(ppi)=mer(ppi)
   	  prevtally(ppi)=tally(ppi)
   	 enddo
   	 prevfrag=frag
   	 endif

!-------VI.L.1. CLUSTER SNAPSHOTS - MATCHING TARGET TIME AND MER
         
  	 if (RPS .eq. 1) then
 	 qik = ctme - (1d-5)
 	 qok = ctme + (1d-5)
 	        
    	 sphf = '(a11,f10.5,a11,i7,a20,i8,a8,i8)'
   	 sf   = '(i5,a3,a7,i5,3f8.3)'
   	 sg   = '(i5,a4,a6,i5,3f8.3)'
 	 sh   = '(3f13.5)'

 	 do np = 1,5                         ! roll call of target time.
 	  if ( ( qik .lt. gt( np ) )  .and.
     &         ( qok .gt. gt( np ) ) ) then
 	   do nth = 1,le( np )               ! roll call of target mers.
 	    nin = gmin( np,nth )             ! lower limit of chosen range.
 	    nfn = gmax( np,nth )             ! upper limit of chosen range.
 	    nap = snap( np,nth )             ! lower limit of chosen range.
 	
    	    do ffi = 1,frag
    	     FRGPOOL( ffi ) = .TRUE.
    	    enddo
   	
   	    LABEL = .FALSE.
   	    SNAPD = .FALSE.

!-------VI.L.2. CLUSTER SNAPSHOTS - WRITE TOTAL ATOMS

   	    do mmi = 1,pop           
 	     if ( ( tally( mmi ) .ne. 0 ) .and. ( mmi .eq. nap ) ) then 
    	      do ppi = 1,pop  
 	       if ( mer( ppi ) .eq. mmi ) then
 	        do ffj = 1,frag  
 	         if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &                ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	      
 	          FRGPOOL( ffj ) = .FALSE.  
	          typcnt1 = 0              
	          typcnt2 = 0             
	          typcnt3 = 0             
	          typcnt4 = 0               
	          typcnt5 = 0               
	          typcnt6 = 0               
	          typcnt8 = 0               

             
	          do ppj = 1,pop                                               
    	           if ( nuc( ppj ) .eq. ffj ) then
	            if     ( molty( ppj ) .eq. 1 ) then
	             typcnt1 = typcnt1 + 1
 	             l = 4
	            elseif ( molty( ppj ) .eq. 2 ) then
	             typcnt2 = typcnt2 + 1
 	             l = 9
	            elseif ( molty( ppj ) .eq. 3 ) then
	             typcnt3 = typcnt3 + 1
 	             l = 6
	            elseif ( molty( ppj ) .eq. 4 ) then
	             typcnt4 = typcnt4 + 1
 	             l = 5
	            elseif ( molty( ppj ) .eq. 5 ) then
	             typcnt5 = typcnt5 + 1
 	             l = 3
	            elseif ( molty( ppj ) .eq. 6 ) then
	             typcnt6 = typcnt6 + 1
 	             l = 5
	            elseif ( molty( ppj ) .eq. 8 ) then
	             typcnt8 = typcnt8 + 1
 	             l = 10
  	            endif
 	           endif
 	          enddo
 	          skips=typcnt1*4+typcnt2*9+typcnt3*6+
     &                  typcnt4*5+typcnt5*3+typcnt6*5+
     &                  typcnt8*10
 	          rmd = 5000
 	          mxv = np*10
 	          rmc = rmd + mxv + nth
 	          write(rmc,sphf) "TIME (ns) =",ctme,"MER = ",mmi,
     &                         "MER RANGE =",nin,"to",nfn
 	          write(rmc,*) skips
 	          LABEL = .TRUE.
 	          EXIT
 	         endif
 	        enddo
 	       endif 
 	
 	       if ( LABEL .eqv. .TRUE. ) then
 	        EXIT
 	       endif
 	      enddo
 	     endif
 	
 	     if ( LABEL .eqv. .TRUE. ) then
 	      EXIT
 	     endif
 	    enddo

    	    do ffi = 1,frag
    	     FRGPOOL( ffi ) = .TRUE.
    	    enddo
                 
!-------VI.L.3. CLUSTER SNAPSHOTS - ITERATE ALL MER SIZE

   	    do mmi  = 1,pop           
 	     if ( ( tally( mmi ) .ne. 0 ) .and. ( mmi .eq. nap ) ) then 

!-------VI.L.4. CLUSTER SNAPSHOTS - SEARCHING FRAGS INVOLVED IN CURRENT CLUSTER SIZE

    	      do ppi = 1,pop  
 	       if ( mer( ppi ) .eq. mmi ) then
 	        do ffj = 1,frag  
 	         if ( ( nuc( ppi ) .eq. ffj ) .and.                
     &                ( FRGPOOL( ffj ) .eqv. .TRUE. ) ) then
 	      
 	          FRGPOOL( ffj ) = .FALSE.  
	          typcnt1 = 0              
	          typcnt2 = 0             
	          typcnt3 = 0             
                  typcnt4 = 0               
	          typcnt5 = 0               
	          typcnt6 = 0               
	          typcnt8 = 0               
                  cq = 0
                  cp = 0
                          
!-------VI.L.5. CLUSTER SNAPSHOTS - INDEXING MOLECULES INVOLVED IN LOCATED FRAG
                  
	          do ppj = 1,pop                                               
    	           if ( nuc( ppj ) .eq. ffj ) then
	            if     ( molty( ppj ) .eq. 1 ) then
	             cq = cq + 1
	             typcnt1 = typcnt1 + 1
	             l = 4
	            elseif ( molty( ppj ) .eq. 2 ) then
	             cq = cq + 1
	             typcnt2 = typcnt2 + 1
	             l = 9
	            elseif ( molty( ppj ) .eq. 3 ) then
	             cq = cq + 1
	             typcnt3 = typcnt3 + 1
	             l = 6
	            elseif ( molty( ppj ) .eq. 4 ) then
	             cq = cq + 1
	             typcnt4 = typcnt4 + 1
	             l = 5
	            elseif ( molty( ppj ) .eq. 5 ) then
	             cq = cq + 1
	             typcnt5 = typcnt5 + 1
	             l = 3
	            elseif ( molty( ppj ) .eq. 6 ) then
	             cq = cq + 1
	             typcnt6 = typcnt6 + 1
 	             l = 5
	            elseif ( molty( ppj ) .eq. 8 ) then
	             cq = cq + 1
	             typcnt8 = typcnt8 + 1
 	             l = 10
  	            endif

!-------VI.L.6. CLUSTER SNAPSHOTS - COMMENCE WRITING OF STRUCTURE FILE

   	  k = molty(ppj)
   	 if ( k .eq. 1 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"SOL"," OW",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sf) cq,"SOL","HW1",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sf) cq,"SOL","HW2",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1
   	 write(rmc,sf) cq,"SOL"," MW",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	
   	 elseif ( k .eq. 2 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09A",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09B",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09C",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09D",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09E",cp,x(ppj,k,5),y(ppj,k,5),z(ppj,k,5)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09F",cp,x(ppj,k,6),y(ppj,k,6),z(ppj,k,6)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09G",cp,x(ppj,k,6),y(ppj,k,6),z(ppj,k,6)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09H",cp,x(ppj,k,7),y(ppj,k,7),z(ppj,k,7)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NON","C09I",cp,x(ppj,k,8),y(ppj,k,8),z(ppj,k,8)

   	 elseif ( k .eq. 3 ) then
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","CAA",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","CAB",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","CAC",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","CAD",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","OAE",cp,x(ppj,k,5),y(ppj,k,5),z(ppj,k,5)
   	 cp=cp+1
   	 write(rmc,sg) cq,"BUTA","HAJ",cp,x(ppj,k,6),y(ppj,k,6),z(ppj,k,6)

   	 elseif ( k .eq. 4 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"NH3","N1",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NH3","H1",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NH3","H2",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NH3","H3",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	 cp=cp+1
   	 write(rmc,sf) cq,"NH3","DN",cp,x(ppj,k,5),y(ppj,k,5),z(ppj,k,5)

   	 elseif ( k .eq. 5 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"MET","CBA",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sf) cq,"MET","OBB",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sf) cq,"MET","HBC",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
          
   	 elseif ( k .eq. 6 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"ACT","CQA",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1
   	 write(rmc,sf) cq,"ACT","CQB",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1
   	 write(rmc,sf) cq,"ACT","OQA",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1
   	 write(rmc,sf) cq,"ACT","OQB",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	 cp=cp+1
   	 write(rmc,sf) cq,"ACT","HQA",cp,x(ppj,k,5),y(ppj,k,5),z(ppj,k,5)

   	 elseif ( k .eq. 8 ) then
   	 cp=cp+1
   	 write(rmc,sf) cq,"OCTA","CAA",cp,x(ppj,k,1),y(ppj,k,1),z(ppj,k,1)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAB",cp,x(ppj,k,2),y(ppj,k,2),z(ppj,k,2)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAC",cp,x(ppj,k,3),y(ppj,k,3),z(ppj,k,3)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAD",cp,x(ppj,k,4),y(ppj,k,4),z(ppj,k,4)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAE",cp,x(ppj,k,5),y(ppj,k,5),z(ppj,k,5)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAF",cp,x(ppj,k,6),y(ppj,k,6),z(ppj,k,6)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAG",cp,x(ppj,k,7),y(ppj,k,7),z(ppj,k,7)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","CAH",cp,x(ppj,k,8),y(ppj,k,8),z(ppj,k,8)
   	 cp=cp+1                                                         
   	 write(rmc,sf) cq,"OCTA","OAE",cp,x(ppj,k,9),y(ppj,k,9),z(ppj,k,9)
   	 cp=cp+1
   	 write(rmc,sf) cq,"OCTA","HAJ",
     &                 cp,x(ppj,k,10),y(ppj,k,10),z(ppj,k,10)

   	 endif
   	           endif
   	          enddo
   	          SNAPD = .TRUE.
   	
   	          if ( SNAPD .eqv. .TRUE. ) then
   	           EXIT
   	          endif
   	
   	         endif                       ! only molec contained in current frag/cluster may proceed.
   	        enddo                        ! end of searching frags with current cluster size.
                write(rmc,sh) sybx,sybx,sybx 
   	       endif                         ! only molec contained in current mer size may proceed.
   	
   	       if ( SNAPD .eqv. .TRUE. ) then
   	        EXIT
   	       endif
   	      enddo                          ! end of iterating through all molecs. 
   	     endif
   	
   	     if ( SNAPD .eqv. .TRUE. ) then
   	      EXIT
   	     endif
   	    enddo                            ! end of iterating all cluster sizes.
   	   enddo                             ! roll call of target cluster sizes.
   	
   	  endif                              ! only time matching with target time may proceed. 
  	 enddo                               ! roll call of target time.
  	 endif

!-------VI.M. MOLECULE REPORT - WRITE MOLECULES INFORMATION, input in growth.f

c	 do cc=1,pop
c	  write(8000,*) "mol",fr,cc,nuc(cc),mer(cc)
c	 enddo

!-------VI.N. NUCLEI REPORT - WRITE NUCLEI INFORMATION

   	 if (RPTR .eq. 1) then
	 ggg=9000+fr
 	 write(ggg,'(a6,i8,a15,f10.4)') "FRAME",fr,"TIME (ns)",ctme
   	 write(ggg,*)
   	 write(ggg,*) "AUDITING"
   	 audsum=0
	 audsum2=0
   	 do mercnt=1,pop
 	  if ( mercnt .eq. 1 ) then
c 	   write(8,*) ".",fr,tally(mercnt)
 	  endif
  	  if ( tally(mercnt) .ne. 0 ) then
  	   write(ggg,*) mercnt,tally(mercnt),tally(mercnt)/mercnt
   	   audsum=audsum+tally(mercnt)
   	   audsum2=audsum2+INT(tally(mercnt)/mercnt)
   	  endif
  	 enddo

   	 write(ggg,*) 
  	 write(ggg,*) "TOTAL MOLECULES:",audsum
  	 write(ggg,*) "TOTAL TALLY:",audsum2
   	 write(ggg,*) "TOTAL FRAGMENTS:",frag
   	 write(ggg,*) "RECOVERED:",recov
   	
   	 do i=1,frag
   	  typcnt1=0
   	  typcnt2=0
   	  typcnt3=0
   	  typcnt4=0
   	  typcnt5=0
   	  typcnt6=0
   	  typcnt8=0
   	  write(ggg,*) 
   	  write(ggg,*) "****NUCLEUS",i
   	  do cc=1,pop
   	   if ( nuc(cc) .eq. i ) then
   	    k=molty(cc)
   	   if ( k .eq. 1 ) then
   	    typcnt1=typcnt1+1
   	   elseif ( k .eq. 2 ) then
   	    typcnt2=typcnt2+1
   	   elseif ( k .eq. 3 ) then
   	    typcnt3=typcnt3+1
   	   elseif ( k .eq. 4 ) then
   	    typcnt4=typcnt4+1
   	   elseif ( k .eq. 5 ) then
   	    typcnt5=typcnt5+1
   	   elseif ( k .eq. 6 ) then
   	    typcnt6=typcnt6+1
   	   elseif ( k .eq. 8 ) then
   	    typcnt8=typcnt8+1
   	   endif
   	  endif
   	  enddo
   	 clustersize=typcnt1+typcnt2+typcnt3+typcnt4+typcnt5+
     &               typcnt6+typcnt8
   	 skips=typcnt1*4+typcnt2*9+typcnt3*6+typcnt4*5+
     &         typcnt5*3+typcnt6*5+typcnt8*10
   	 write(ggg,*) "Size/WNBAMC/atoms",clustersize,"|",
     &    typcnt1,typcnt2,typcnt3,typcnt4,typcnt5,
     &    typcnt6,typcnt8,"|",skips
   	 do cc=1,pop
   	  if ( nuc(cc) .eq. i ) then
   	   k=molty(cc)
   	   if ( k .eq. 1 ) then
   	    write(ggg,*) "O",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "H",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "H",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "M",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	   elseif ( k .eq. 2 ) then
   	    write(ggg,*) "C",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "C",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "C",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "C",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	    write(ggg,*) "C",x(cc,k,5),y(cc,k,5),z(cc,k,5)
   	    write(ggg,*) "C",x(cc,k,6),y(cc,k,6),z(cc,k,6)
   	    write(ggg,*) "C",x(cc,k,7),y(cc,k,7),z(cc,k,7)
   	    write(ggg,*) "C",x(cc,k,8),y(cc,k,8),z(cc,k,8)
   	    write(ggg,*) "C",x(cc,k,9),y(cc,k,9),z(cc,k,9)
   	   elseif ( k .eq. 3 ) then
   	    write(ggg,*) "C",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "C",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "C",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "C",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	    write(ggg,*) "O",x(cc,k,5),y(cc,k,5),z(cc,k,5)
   	    write(ggg,*) "H",x(cc,k,6),y(cc,k,6),z(cc,k,6)
   	   elseif ( k .eq. 4 ) then
   	    write(ggg,*) "N",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "H",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "H",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "H",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	    write(ggg,*) "M",x(cc,k,5),y(cc,k,5),z(cc,k,5)
   	   elseif ( k .eq. 5 ) then
   	    write(ggg,*) "C",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "O",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "H",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	   elseif ( k .eq. 6 ) then
   	    write(ggg,*) "C",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "C",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "O",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "O",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	    write(ggg,*) "H",x(cc,k,5),y(cc,k,5),z(cc,k,5)
   	   elseif ( k .eq. 8 ) then
   	    write(ggg,*) "C",x(cc,k,1),y(cc,k,1),z(cc,k,1)
   	    write(ggg,*) "C",x(cc,k,2),y(cc,k,2),z(cc,k,2)
   	    write(ggg,*) "C",x(cc,k,3),y(cc,k,3),z(cc,k,3)
   	    write(ggg,*) "C",x(cc,k,4),y(cc,k,4),z(cc,k,4)
   	    write(ggg,*) "C",x(cc,k,5),y(cc,k,5),z(cc,k,5)
   	    write(ggg,*) "C",x(cc,k,6),y(cc,k,6),z(cc,k,6)
   	    write(ggg,*) "C",x(cc,k,7),y(cc,k,7),z(cc,k,7)
   	    write(ggg,*) "C",x(cc,k,8),y(cc,k,8),z(cc,k,8)
   	    write(ggg,*) "O",x(cc,k,9),y(cc,k,9),z(cc,k,9)
   	    write(ggg,*) "H",x(cc,k,10),y(cc,k,10),z(cc,k,10)
   	   endif
   	  endif
   	 enddo
   	 enddo
   	 endif

!-------VI.O. MER SIZE REPORT - SUMMARY OF MOLECULES FOR EACH CLUSTER SIZE

c	 do i=1,pop
c         k=molty(i)
c         fff=10000+mer(i)
c 	  if ( k .eq. 1 ) then
c  	   write(fff,*) "O",x(i,k,1),y(i,k,1),z(i,k,1)
c   	   write(fff,*) "H",x(i,k,2),y(i,k,2),z(i,k,2)     
c   	   write(fff,*) "H",x(i,k,3),y(i,k,3),z(i,k,3)
c   	   write(fff,*) "M",x(i,k,4),y(i,k,4),z(i,k,4)
c  	  elseif ( k .eq. 2 ) then
c  	   write(fff,*) "C",x(i,k,1),y(i,k,1),z(i,k,1)
c  	   write(fff,*) "C",x(i,k,2),y(i,k,2),z(i,k,2)
c   	   write(fff,*) "C",x(i,k,3),y(i,k,3),z(i,k,3)
c   	   write(fff,*) "C",x(i,k,4),y(i,k,4),z(i,k,4)
c   	   write(fff,*) "C",x(i,k,5),y(i,k,5),z(i,k,5)
c   	   write(fff,*) "C",x(i,k,6),y(i,k,6),z(i,k,6)
c   	   write(fff,*) "C",x(i,k,7),y(i,k,7),z(i,k,7)
c   	   write(fff,*) "C",x(i,k,8),y(i,k,8),z(i,k,8)
c   	   write(fff,*) "C",x(i,k,9),y(i,k,9),z(i,k,9)
c  	  elseif ( k .eq. 3 ) then
c   	   write(fff,*) "C",x(i,k,1),y(i,k,1),z(i,k,1)
c   	   write(fff,*) "C",x(i,k,2),y(i,k,2),z(i,k,2)
c   	   write(fff,*) "C",x(i,k,3),y(i,k,3),z(i,k,3)
c   	   write(fff,*) "C",x(i,k,4),y(i,k,4),z(i,k,4)
c   	   write(fff,*) "O",x(i,k,5),y(i,k,5),z(i,k,5)
c   	   write(fff,*) "H",x(i,k,6),y(i,k,6),z(i,k,6)
c         elseif ( k .eq. 4 ) then
c          write(fff,*) "N",x(i,k,1),y(i,k,1),z(i,k,1)
c 	   write(fff,*) "H",x(i,k,2),y(i,k,2),z(i,k,2)
c 	   write(fff,*) "H",x(i,k,3),y(i,k,3),z(i,k,3)
c 	   write(fff,*) "H",x(i,k,4),y(i,k,4),z(i,k,4)
c  	   write(fff,*) "M",x(i,k,5),y(i,k,5),z(i,k,5)
c         elseif ( k .eq. 5 ) then
c          write(fff,*) "AR",x(i,k,1),y(i,k,1),z(i,k,1)
c  	  endif
c        enddo

	enddo          ! end frame scans

! ########################################################################################################

	call system_clock (t2, clock_rate, clock_max)
	SECONDS=real(t2-t1)/real(clock_rate)
	write(*,*) " "
	write(*,'(a21,i3,a9,f7.3,a10)')'      (TIME ELAPSED:',
     &                          INT(SECONDS/60),'MINUTE/S',
     &                          MOD(SECONDS,60d0),'SECOND/S)'
	write(*,*) " "

! ########################################################################################################

! SUBROUTINE 

	CONTAINS
	 subroutine UNDO_TRANSLATE(pi)
	  integer pi

	  pool(pi)=.TRUE.
	  RNDTAG(pi)=.FALSE.
	  if ( molty(pi) .eq. 1 ) then
	   lp=4
	  elseif ( molty(pi) .eq. 2 ) then
	   lp=9
	  elseif ( molty(pi) .eq. 3 ) then
	   lp=6
	  elseif ( molty(pi) .eq. 4 ) then
	   lp=5
	  elseif ( molty(pi) .eq. 5 ) then
	   lp=3
	  elseif ( molty(pi) .eq. 6 ) then
	   lp=5
	  elseif ( molty(pi) .eq. 7 ) then
	   lp=1
	  elseif ( molty(pi) .eq. 8 ) then
	   lp=1
	  endif

	  if ( pxmore(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    x(pi,molty(pi),pb)=x(pi,molty(pi),pb)+sybx
	   enddo
	  endif
	  if ( pymore(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    y(pi,molty(pi),pb)=y(pi,molty(pi),pb)+sybx
	   enddo
	  endif
	  if ( pzmore(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    z(pi,molty(pi),pb)=z(pi,molty(pi),pb)+sybx
	   enddo
	  endif
	  if ( pxless(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    x(pi,molty(pi),pb)=x(pi,molty(pi),pb)-sybx
	  enddo
	  endif
	  if ( pyless(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    y(pi,molty(pi),pb)=y(pi,molty(pi),pb)-sybx
	   enddo
	  endif
	  if ( pzless(pi) .eqv. .TRUE. ) then
	   do pb=1,lp
	    z(pi,molty(pi),pb)=z(pi,molty(pi),pb)-sybx
	   enddo
	  endif
	 end subroutine UNDO_TRANSLATE

! ########################################################################################################

	end
