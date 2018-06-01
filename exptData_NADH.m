function [timeROF, timeOCR, timeIncEnDem, timeFR, NADH_Glc_ROF_FC,...
    NADH_Glc_OF_FC, NADH_Glc_IncEnDem_FC, NADH_Glc_FR_FC,...
    NADH_Glc_ROF_cal, NADH_Glc_OF_cal, NADH_Glc_IncEnDem_cal,...
    NADH_Glc_FR_cal] = exptData_NADH 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental data from Padova WT cells (median)
% Taken from "Autofluo Padova - Niamh.xls" - median of WT
% NADH_Glc_ROF_FC = NADH (FC) in response to Rotenone and Oligo 
% NADH_Glc_OF_FC = NADH (FC) in response to Oligo and FCCP
% NADH_Glc_FR_FC = NADH (FC) in response to FCCP and Rotenone
% NADH_Glc_IncEnDem_FC = NADH in response to Gramicidin 
% NADH_Glc_ROF_cal = calibrated NADH(FC) ROF
% NADH_Glc_OF_cal = calibrated NADH(FC) OF
% NADH_Glc_FR_cal = calibrated NADH(FC) FR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% False time to align with simulations
timeROF = 5:1/12:20; % Rotenone expt
timeOCR = [13:1/12:21.2 34.25:1/12:41]; % Seahorse expt (Oligo + FCCP)
timeIncEnDem = -2:1/12:23;
timeFR = [5:1/12:12.2 37.2:1/12:45];

% Import FC values
NADH_Glc_ROF_FC = [1.029659574
1.022481829
1.010089784
1.01520858
1.010851651
1.014619829
1.012022445
1.016782756
1.00650344
1.009878258
0.997901684
1.010633466
0.99636096
1.005914615
0.997346804
1.001117262
1.008518325
0.995018022
1.006895262
0.992601407
0.994937942
0.997042321
1.001885934
0.997143795
0.991357663
1.001052426
1.003082034
0.989446136
0.997800927
0.994946379
1.004595553
0.992415356
0.993315282
1.003142642
0.9859601
1.00043204
1.001676301
0.994874764
0.995893417
0.996840005
1.008190766
0.991511071
0.98045105
0.994635242
1.000733004
0.994692955
1.000181915
0.993879147
0.989124402
0.993595277
0.987387064
1.00548718
0.998708955
0.99735792
0.991694354
0.983825622
0.985121287
0.995867446
0.991688987
0.996989236
1.001262824
1.002248818
1.006575894
1.050212512
1.114192587
1.185830147
1.265501881
1.290455579
1.320575876
1.352243984
1.365542259
1.363877186
1.365422825
1.384010269
1.388678256
1.39456559
1.386231363
1.390465946
1.412158474
1.401418013
1.438394203
1.415623768
1.414470026
1.418717679
1.427896259
1.431999993
1.430884341
1.440939134
1.44326858
1.455644598
1.43895216
1.461136288
1.455888504
1.446683187
1.477534501
1.477679788
1.480845825
1.475101751
1.479173807
1.491931267
1.478002058
1.485597632
1.478470069
1.490025446
1.485814325
1.497216567
1.484779939
1.501529528
1.512263224
1.501560505
1.496387492
1.50582189
1.512878879
1.509573375
1.496984299
1.511042543
1.49598128
1.505804989
1.503369948
1.501348853
1.550039595
1.549447764
1.500542081
1.520339266
1.536451532
1.51377425
1.519027376
1.514195994
1.520839191
1.523824705
1.51047817
1.519573497
1.513355563
1.539149872
1.527262371
1.521900907
1.530829139
1.560527133
1.507904304
1.526554069
1.555867946
1.528777957
1.518290125
1.535335477
1.524679857
1.533556877
1.56153966
1.546525611
1.556543248
1.52103081
1.536354137
1.549319503
1.546048448
1.558091693
1.540960111
1.541195794
1.531304526
1.553041996
1.537757599
1.523405456
1.545106419
1.562405883
1.550504976
1.525893719
1.539858537
1.557268604
1.542752396
1.546940828
1.555577388
1.521997224
1.55646316
1.562544133
1.539198161
1.534793845
1.542573585
1.549717714
1.563640655
1.572740689
1.545979479
1.556517985
1.534631302
];

NADH_Glc_OF_FC = [1.037871128
1.035320671
1.031075289
1.039898022
1.018841604
1.021037375
1.022468864
1.017794935
1.028897482
1.006717855
1.00318879
1.010324445
1.018573577
1.00160518
1.004678852
0.996715475
1.003341802
1.004763653
1.008884858
0.993699833
0.994096857
0.999138116
0.993733755
1.004480644
1.01196712
0.988085415
0.989762329
0.992604738
0.98676167
0.989916316
0.998171984
0.990147096
0.989570702
1.002783126
1.002786517
0.997475865
1.002064249
1.000889781
0.999693919
0.998168352
0.998483425
0.968668032
0.999550556
0.984341779
0.99284402
0.985223611
0.988510262
0.996005524
0.980411324
0.992915916
1.001143718
0.992350403
0.990110745
0.980508694
0.984170288
0.987609594
0.986880685
0.994852401
0.992606203
0.995324369
0.983747257
0.991225081
0.982792687
0.987325699
1.001795871
1.023476516
1.056121548
1.111025487
1.20357131
1.263282425
1.326371561
1.367594846
1.379210288
1.414326931
1.426038561
1.413360347
1.448731047
1.467199509
1.465385578
1.478820788
1.472592751
1.509094815
1.478033668
1.469836299
1.482332791
1.462424027
1.471834352
1.45410358
1.455728901
1.46442439
1.474302707
1.478054405
1.515411465
1.487301862
1.473123028
1.463615524
1.472317489
1.48318025
1.461882489
1.47068338
1.472716994
1.491019427
1.485552098
1.462008084
1.482015782
1.492274618
1.496319189
1.483915725
1.496802481
1.459074292
1.481604735
1.465147192
1.488624562
1.499218787
1.494045593
1.47470597
1.470858815
1.480616223
1.496802481
1.474411703
1.403322541
0.930255982
0.755284946
0.729051742
0.690204587
0.68343129
0.680561902
0.658687951
0.666714689
0.673545345
0.654458528
0.660608852
0.660009032
0.658865345
0.641269024
0.676892925
0.668674107
0.675681546
0.65820321
0.67875237
0.666146569
0.669865563
0.650023654
0.669353507
0.680667389
0.694043733
0.667647291
0.655429604
0.674310158
0.637642405
0.68878135
0.640788093
0.673766254
0.651603225
0.66243389
0.686863414
0.683145131
0.717681098
0.643633364
0.691019205
0.663870505
0.684183383
0.678051334
0.6500359
0.661698971
0.684229448
0.659063661
0.648234951
0.660229049
0.666108653
0.66531705
0.661446304
0.647489017
0.645866621
0.669969044
0.630341128
0.629806844
0.619842869
0.640160527
0.639037923
0.635444516
    ];

NADH_Glc_FR_FC = [1.050086133
1.054765248
1.023376991
1.031519437
1.030454548
1.023355105
1.011760801
0.997922888
1.013442729
1.014905394
1.00785503
1.011318741
1.004895709
1.003033126
1.019864767
0.998084366
1.009067429
1.001034413
1.008704242
1.002099136
0.986508751
1.012296413
0.994895335
0.992121312
1.005300226
1.009673838
0.994626365
0.983918296
1.005627304
1.000835079
1.01460643
0.990118
0.991967255
0.986173994
0.989034582
0.993528669
0.992092719
0.994804781
0.986416872
0.987605248
0.973451648
0.979076134
0.977803502
0.975503928
0.990800955
1.003188125
0.972180968
0.979503786
0.994049216
1.007792013
0.993928909
0.994824428
0.979847837
0.981472855
0.995930851
0.986146282
0.993789373
0.988241527
1.009893689
0.980297757
0.985204954
0.93108252
0.874880882
0.848690732
0.776286573
0.746467361
0.742061283
0.737042618
0.704181609
0.712607226
0.707671161
0.702198223
0.673252285
0.656585757
0.673918051
0.681891355
0.703699794
0.70405663
0.692749323
0.691790429
0.670660848
0.661684707
0.670336503
0.667880903
0.688677584
0.656868426
0.669337833
0.682539735
0.68251234
0.687165493
0.674927276
0.665017697
0.672527725
0.685928559
0.67354581
0.673022194
0.670372069
0.669968398
0.705246931
0.676231928
0.68280518
0.68307644
0.667174719
0.667702936
0.65446703
0.67260378
0.663785754
0.673371673
0.692599715
0.671876042
0.659701032
0.673501942
0.638719075
0.639394826
0.63696124
0.648596424
0.669763409
0.659227656
0.661579801
0.653655747
0.64506385
0.663820406
0.668595129
0.700941665
0.74555709
0.835413519
0.929374182
1.010144887
1.001515194
1.030857727
1.088585217
1.101930933
1.100261727
1.136410127
1.133735895
1.129390016
1.137264898
1.118393429
1.105920069
1.089432071
1.114000414
1.122448922
1.115697349
1.092754837
1.10371319
1.087216328
1.106594135
1.082035744
1.077681587
1.103921606
1.081016522
1.065140817
1.075707184
1.073785071
1.069663233
1.059323864
1.075522475
1.065693651
1.063761161
1.094659853
1.056796306
1.056654033
1.057628904
1.051112639
1.047408536
1.065066524
1.029139661
1.019250441
1.032109893
1.04093055
1.04741014
1.02217911
1.022787368
1.03190793
1.038429613
1.008337684
1.062505092
1.018841923
1.041685617
1.035955923
1.052140006
];

NADH_Glc_IncEnDem_FC = [1.052269029
1.045595557
1.02618614
1.042071079
1.025684112
1.021675112
1.019971843
1.005613386
1.003064474
1.007592057
1.010478091
1.019599611
1.008172518
1.014826586
1.002317379
1.009836151
0.99704107
1.003010575
1.008665024
0.990254922
1.003343906
0.997832522
0.996633127
0.999631783
0.986108448
0.988776258
0.993997404
0.994505208
0.98922281
0.997405719
0.988732959
0.997046944
1.007081371
0.999200327
0.988301043
0.986020911
0.986322895
0.98848506
0.991829972
0.998623196
0.984639325
0.992309793
0.992278604
0.974735031
0.978452612
0.990139304
0.991277246
0.989733512
1.007772682
0.981695048
0.999196696
0.984908826
0.984834897
0.99135553
0.995465715
0.99105819
0.986501089
0.987806166
0.979014792
1.000105013
0.9893158
0.987055749
0.964658953
0.954077549
0.94821652
0.956138755
0.944615591
0.932449968
0.931188524
0.947946836
0.920975978
0.92693666
0.919488083
0.896338774
0.885749314
0.879535813
0.897469334
0.877068262
0.86317489
0.866050309
0.841120374
0.832720858
0.862877397
0.83316166
0.837426585
0.826078436
0.848732343
0.832578351
0.832549879
0.811481299
0.802030369
0.833905522
0.799559413
0.804193498
0.816215795
0.802487134
0.786355435
0.778343363
0.791888851
0.798524135
0.784534698
0.773147004
0.788461698
0.771957526
0.787519493
0.774867125
0.751028289
0.756635617
0.773986927
0.738671816
0.766163222
0.736765771
0.765527282
0.743428989
0.727271644
0.770314091
0.74988516
0.728497589
0.716118845
0.707879615
0.720049799
0.7314932
0.71381911
0.70224986
0.714405748
0.728930746
0.718380228
0.6923447
0.7187507
0.706383735
0.689649302
0.7073328
0.690819995
0.707513164
0.696399209
0.68495009
0.702678127
0.699017936
0.700831143
0.691583068
0.687732009
0.689991413
0.69040819
0.677045436
0.71140849
0.704582795
0.667698385
0.675334881
0.680265585
0.680588791
0.689634467
0.675986397
0.673027541
0.657046783
0.643348756
0.665016023
0.646201839
0.67394173
0.655895199
0.663884535
0.656035364
0.661325041
0.66725749
0.663142906
0.668036521
0.651838062
0.662558876
0.661151846
0.641204405
0.644493341
0.647570818
0.673202938
0.638277503
0.651468935
0.633740643
0.645576823
0.632744026
0.654764449
0.642456882
0.66979029
0.649944422
0.637739411
0.654904388
0.647789408
0.613943748
0.633615686
0.617047964
0.61662121
0.64133322
0.636616963
0.631194341
0.635192447
0.618975452
0.628372487
0.615808684
0.628004358
0.636426356
0.612426791
0.615998095
0.615004539
0.62042178
0.611677923
0.604221458
0.603842618
0.617532845
0.600945605
0.60912093
0.625752973
0.632810623
0.606024505
0.616842081
0.615968625
0.580747732
0.606118234
0.606783027
0.599861348
0.632417225
0.60998816
0.612345515
0.612053395
0.589029237
0.5955058
0.592639583
0.612997296
0.62064105
0.609980715
0.578278249
0.58269531
0.578642905
0.567990166
0.618035935
0.62125143
0.57001961
0.570625518
0.572634964
0.602497189
0.590651555
0.585846321
0.613728509
0.58958312
0.574665646
0.596043356
0.565408574
0.586728685
0.57837122
0.585364174
0.55524724
0.554994364
0.543270564
0.578044652
0.59274912
0.611277146
0.574843242
0.60570549
0.570040642
0.563840162
0.565446279
0.548864195
0.583520362
0.552300352
0.530350472
0.601396904
0.542907499
0.553265973
0.534849201
0.548323432
0.558546918
0.55524724
0.545964661
0.560305347
0.54398096
0.536641032
0.524926854
0.551081492
0.517680956
0.557382209
0.541700654
0.558508027
0.53861362
0.52090134
0.55914689
0.55376851
0.530715478
0.493153146
0.535375661
0.525838511
0.522846628
0.527878906
0.512073234
0.527902146
0.55889288
0.522981558
0.525585428
0.535892158
0.552847464
0.543380833
0.533196179
0.551255546
0.510587456
0.501762855
0.539787962
];

% Import calibrated values
% (pasted into MATLAB first to list them as a row instead of a column)
NADH_Glc_ROF_cal = [1.04017608900000;1.03489243700000;1.01485311100000;1.01911926400000;1.01860663600000;1.02480774400000;1.01474045200000;1.02433882000000;1.01049631000000;1.01453522000000;0.997302014000000;1.01382585400000;0.995273986000000;1.00906956800000;0.996324210000000;1.00130952200000;1.01332285400000;0.993015487000000;1.00974965600000;0.991140038000000;0.993294331000000;0.995432973000000;1.00285270400000;0.995833391000000;0.986801987000000;1.00195677200000;1.00425106100000;0.985006586000000;0.997090359000000;0.991315600000000;1.00791773700000;0.988727849000000;0.989368833000000;1.00441692600000;0.980905887000000;1.00047318700000;1.00270455500000;0.992863258000000;0.993190638000000;0.995730636000000;1.01095451100000;0.987730904000000;0.975254139000000;0.991122927000000;1.00097464800000;0.992213514000000;1.00026928200000;0.991337224000000;0.985880599000000;0.989988973000000;0.980955144000000;1.00904888100000;0.998201384000000;0.996439794000000;0.988615792000000;0.977154640000000;0.979882211000000;0.992480272000000;0.985996487000000;0.995706608000000;1.00266392400000;1.00378414500000;1.00765511800000;1.07814388400000;1.16401243200000;1.28377929500000;1.37078037100000;1.42981966200000;1.45564949100000;1.46696363000000;1.47825092800000;1.51509482100000;1.51467035800000;1.56343157200000;1.56530943800000;1.57101535900000;1.56297873300000;1.56170640800000;1.58896403700000;1.56629801900000;1.57369371800000;1.61365992700000;1.60458637300000;1.60376122000000;1.64236122000000;1.60835626300000;1.63754895100000;1.64281431400000;1.62598137700000;1.65146385200000;1.64506296900000;1.67485369400000;1.68373822400000;1.66303093800000;1.71756895400000;1.71322756100000;1.71082985300000;1.69822594800000;1.70962621600000;1.72736562000000;1.70777786000000;1.75359101500000;1.71363286400000;1.72351727300000;1.73667160200000;1.73854106200000;1.75942644100000;1.73934498900000;1.75605450800000;1.73140995600000;1.72298397000000;1.77228561000000;1.73259733400000;1.74794611900000;1.75544648300000;1.74450608900000;1.78926501500000;1.77483850500000;1.78426583400000;1.80762536000000;1.79547136700000;1.79633667000000;1.80483843800000;1.80441022400000;1.81459522600000;1.81602304900000;1.80848194500000;1.82367368000000;1.79211762900000;1.83285546000000;1.83031326500000;1.81982956500000;1.81833018800000;1.84533981200000;1.82155144900000;1.84997125600000;1.83339621700000;1.87871703800000;1.83416030000000;1.84887039400000;1.86548590000000;1.81645505400000;1.85534805500000;1.81089731700000;1.82017647900000;1.83502052700000;1.87905834900000;1.85919310700000;1.84569623200000;1.82414162700000;1.80696867000000;1.86207580600000;1.82803656900000;1.85121193600000;1.81232470300000;1.84965924800000;1.85764029200000;1.88773115400000;1.83278998200000;1.85321606900000;1.84806970500000;1.88481912100000;1.86458141700000;1.83450332500000;1.83407330600000;1.86038757200000;1.83834031200000;1.84046062500000;1.86434799100000;1.84685258300000;1.88082610700000;1.87620567500000;1.87742640000000;1.82068697600000;1.82727068400000;1.85865228600000;1.82774158300000;1.86936541200000;1.84349459000000;1.84487999200000;1.85161761700000];

NADH_Glc_OF_cal = [1.04847678800000;1.05673341600000;1.04656627600000;1.05569487900000;1.03121046800000;1.02733498500000;1.02987206600000;1.02697763400000;1.04615481700000;1.01129037800000;1.00410597100000;1.02040903600000;1.02829951400000;1.00237972900000;1.00561992500000;0.995463388000000;1.00484929000000;1.01019676900000;1.01108746000000;0.988677280000000;0.991762755000000;0.998809559000000;0.990090408000000;1.00685959000000;1.02494328900000;0.980589306000000;0.986817699000000;0.990560805000000;0.979859984000000;0.985904789000000;0.997748047000000;0.986900691000000;0.981608442000000;1.00501035800000;1.00407965200000;0.996542811000000;1.00313560900000;1.00174872500000;0.999608202000000;0.996030323000000;0.997779619000000;0.956407431000000;0.999242558000000;0.976926596000000;0.986002799000000;0.980249146000000;0.979490385000000;0.993804542000000;0.971320684000000;0.988968616000000;1.00146401100000;0.990576378000000;0.983888880000000;0.973078454000000;0.974217714000000;0.982400622000000;0.981661512000000;0.991149206000000;0.990476109000000;0.992841897000000;0.969030092000000;0.983125127000000;0.975685159000000;0.977706620000000;1.00238758800000;1.04074949500000;1.10750451200000;1.19549263100000;1.40118626600000;1.46509650800000;1.53638836300000;1.64760873700000;1.64232758500000;1.67946480000000;1.69521353300000;1.66481836600000;1.68697951800000;1.72148539900000;1.69125404400000;1.72223976300000;1.75823634100000;1.79346350500000;1.72530787900000;1.70451082600000;1.70843400800000;1.75366540200000;1.73552246300000;1.72343903100000;1.74898584400000;1.69673633900000;1.72343424000000;1.73954595900000;1.75389609800000;1.79993246800000;1.79027104300000;1.73733945300000;1.75285769600000;1.73263561200000;1.73236315600000;1.73565230200000;1.71754739200000;1.78030140200000;1.75754547600000;1.74850182200000;1.77757572900000;1.73305935800000;1.75371595800000;1.76810069900000;1.74698519100000;1.73907825300000;1.75127459000000;1.68399416800000;1.76027015200000;1.76439627400000;1.74280732600000;1.70184247300000;1.78043490800000;1.73967579800000;1.75775799000000;1.74028058100000;1.51611735400000;0.877324212000000;0.551535401000000;0.547998026000000;0.500793304000000;0.428482890000000;0.452569408000000;0.433452091000000;0.450891536000000;0.451040079000000;0.462282645000000;0.455479485000000;0.440197627000000;0.445249174000000;0.446781503000000;0.464676926000000;0.450507453000000;0.466504846000000;0.477607575000000;0.453394190000000;0.437058376000000;0.445062114000000;0.437104694000000;0.447804708000000;0.441847241000000;0.464954116000000;0.455550817000000;0.443287057000000;0.450905080000000;0.439281673000000;0.458146829000000;0.437906003000000;0.449286313000000;0.430124091000000;0.434835144000000;0.465738516000000;0.444506323000000;0.495533731000000;0.421507866000000;0.460614374000000;0.448195449000000;0.459103128000000;0.427643634000000;0.411491577000000;0.415367021000000;0.418599986000000;0.409800867000000;0.387028059000000;0.378535791000000;0.431987787000000;0.399151717000000;0.427677697000000;0.395524298000000;0.399025457000000;0.408379970000000;0.396950770000000;0.383875896000000;0.366914432000000;0.394376791000000;0.382934639000000;0.400856746000000]; 

NADH_Glc_IncEnDem_cal = 0;

NADH_Glc_FR_cal = [[1.07775697800000;1.08719360400000;1.03112279900000;1.05036563700000;1.04370658300000;1.03647247800000;1.01593963400000;0.996515612000000;1.01950117100000;1.02132402400000;1.01195374700000;1.01633041800000;1.00717995300000;1.00515620000000;1.02740420300000;0.997749164000000;1.01307053600000;1.00131297400000;1.01528791400000;1.00274531200000;0.980177432000000;1.01646251100000;0.993284441000000;0.989699335000000;1.00930033900000;1.01418330400000;0.991829123000000;0.972959015000000;1.00779559700000;1.00111553200000;1.02192856200000;0.983865492000000;0.988487157000000;0.980023055000000;0.985110908000000;0.989979565000000;0.989300331000000;0.993650602000000;0.980703042000000;0.981330503000000;0.956821822000000;0.968404834000000;0.967983492000000;0.962477191000000;0.985735522000000;1.00523782300000;0.952792989000000;0.964358320000000;0.992412542000000;1.01088378500000;0.989250536000000;0.991526355000000;0.968607227000000;0.971512675000000;0.993821603000000;0.980796244000000;0.991478542000000;0.984217550000000;1.01427343400000;0.967149252000000;0.978230329000000;0.885022131000000;0.785240447000000;0.736911314000000;0.645714126000000;0.640984036000000;0.621402286000000;0.622691944000000;0.571332005000000;0.564546675000000;0.506082568000000;0.566691038000000;0.525373313000000;0.510887095000000;0.532065915000000;0.517433124000000;0.558498091000000;0.551810054000000;0.542065403000000;0.559527253000000;0.502356702000000;0.503909558000000;0.502367465000000;0.493965980000000;0.509644023000000;0.515980632000000;0.539903461000000;0.505477522000000;0.548681026000000;0.508460557000000;0.501177296000000;0.500443611000000;0.532206269000000;0.503670095000000;0.514068231000000;0.505162185000000;0.541868939000000;0.496306989000000;0.530553589000000;0.483813629000000;0.517433651000000;0.516364656000000;0.510659830000000;0.491297263000000;0.500376710000000;0.537181720000000;0.484715745000000;0.516364962000000;0.539622614000000;0.496706230000000;0.516469479000000;0.500741369000000;0.473977962000000;0.500670185000000;0.463879772000000;0.485081398000000;0.469697121000000;0.470513450000000;0.527630776000000;0.496074879000000;0.508505268000000;0.500331147000000;0.512713547000000;0.549526997000000;0.657185441000000;0.764277756000000;0.903449085000000;1.01505630300000;1.00073712300000;1.04329725600000;1.14803703400000;1.15178035300000;1.16269781900000;1.19953346600000;1.21148041100000;1.18927593200000;1.20641208300000;1.21284345000000;1.14719640000000;1.14000186400000;1.16804147900000;1.18471539000000;1.17703094400000;1.14383579600000;1.15942043600000;1.12652497800000;1.18068725000000;1.12590060700000;1.13494801400000;1.16672715900000;1.12371089300000;1.10059822300000;1.11510626600000;1.12117115000000;1.12449580100000;1.08305115600000;1.11008902000000;1.11334920400000;1.10225493600000;1.14516324400000;1.09336829200000;1.08154892400000;1.07727969300000;1.07474927100000;1.07182686900000;1.09360143800000;1.04013407400000;1.02799231100000;1.04888821100000;1.07007725600000;1.06769287000000;1.03616990300000;1.03335790400000;1.04239253600000;1.05150703200000;1.01104490200000;1.08567901600000;1.02575061700000;1.05638565300000;1.05260483500000;1.07503863000000]];
end