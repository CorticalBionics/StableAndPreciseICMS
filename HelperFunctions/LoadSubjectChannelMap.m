function ChannelMap = LoadSubjectChannelMap(SubjectID)
    valid_subject_ids = {'BCI02', 'BCI03', 'CRS02', 'CRS07', 'CRS08'};
    if strcmpi(SubjectID, 'all')
        ChannelMap = cell(size(valid_subject_ids));
        for s = 1:length(valid_subject_ids)
            ChannelMap{s} = LoadSubjectChannelMap(valid_subject_ids{s});
        end
        ChannelMap = cat(1,ChannelMap{:});
        return
    elseif ~any(strcmpi(valid_subject_ids, SubjectID))
        formatted_string = repmat('\t%s\n', [1, length(valid_subject_ids)]);
        error(['Invalid SubjectID, valid IDs include:\n', formatted_string], valid_subject_ids{:})
    end

    ChannelMap = struct('Subject', [],...
                        'ArrayNames', {}, ...
                        'IsSensory', [], ...
                        'IsMotor', [],...
                        'IsMedial', [], ...
                        'IsLateral', [],...
                        'ArrayLocations', {},...
                        'ArrayOffsets', [],...
                        'ArrayRotations', [],...
                        'ChannelNumbers', {});
    ChannelMap(1).Subject = SubjectID;

    switch SubjectID
        case {'BCI02', 'BCI03', 'CRS07', 'CRS08'} % These all same map but flipped medial/lateral
            if strcmp(SubjectID, 'BCI02')
                ChannelMap(1).ArrayNames = {'MedialMotor', 'MedialSensory', 'LateralMotor', 'LateralSensory'};
                ChannelMap(1).ArrayRotations = [-135, -105, -75, 45];
                ChannelMap(1).IsMedial = logical([1,1,0,0]);
                ChannelMap(1).IsLateral = logical([0,0,1,1]);
            elseif strcmp(SubjectID, 'BCI03')
                ChannelMap(1).ArrayNames = {'MedialMotor', 'MedialSensory', 'LateralMotor', 'LateralSensory'};
                ChannelMap(1).ArrayRotations = [115, 170, -95, -70];
                ChannelMap(1).IsMedial = logical([1,1,0,0]);
                ChannelMap(1).IsLateral = logical([0,0,1,1]);
            elseif strcmp(SubjectID, 'CRS07')
                ChannelMap(1).ArrayNames = {'LateralMotor', 'LateralSensory', 'MedialMotor', 'MedialSensory'};
                ChannelMap(1).ArrayRotations = [70, 75, 10, -10];
                ChannelMap(1).IsMedial = logical([0,0,1,1]);
                ChannelMap(1).IsLateral = logical([1,1,0,0]);
            elseif strcmp(SubjectID, 'CRS08')
                ChannelMap(1).ArrayNames = {'LateralMotor', 'LateralSensory', 'MedialMotor', 'MedialSensory'};
                ChannelMap(1).ArrayRotations = [55, 95, 10, 180];
                ChannelMap(1).IsMedial = logical([0,0,1,1]);
                ChannelMap(1).IsLateral = logical([1,1,0,0]);
            end
            
            ChannelMap(1).IsSensory = logical([0,1,0,1]);
            ChannelMap(1).IsMotor = logical([1,0,1,0]);
            ChannelMap(1).ArrayLocations{1} = [NaN	38	50	59	6	23	22	101	111	NaN;...
                                               33	40	46	64	9	25	24	102	113	128;...
                                               35	43	56	61	17	21	26	103	112	114;...
                                               37	47	55	63	15	14	28	104	115	116;...
                                               39	49	54	2	8	16	30	105	117	118;...
                                               41	48	52	1	11	20	32	106	119	120;...
                                               45	42	58	3	13	27	97	107	121	122;...
                                               34	44	57	4	19	29	99	108	123	124;...
                                               36	51	62	7	10	31	98	109	125	126;...
                                               NaN	53	60	5	12	18	100	110	127	NaN];

            ChannelMap(1).ArrayLocations{2} = [65	NaN	72	NaN	85	91;...
                                               NaN	77	NaN	81	NaN	92;...
                                               67	NaN	74	NaN	87	NaN;...
                                               NaN	79	NaN	82	NaN	93;...
                                               69	NaN	76	NaN	88	NaN;...
                                               NaN	66	NaN	84	NaN	94;...
                                               71	NaN	78	NaN	89	NaN;...
                                               NaN	68	NaN	83	NaN	96;...
                                               73	NaN	80	NaN	90	NaN;...
                                               75	70	NaN	86	NaN	95];

            ChannelMap(1).ArrayLocations{3} = [NaN	166	178	187	134	151	150	229	239	NaN;...
                                               161	168	174	192	137	153	152	230	241	256;...
                                               163	171	184	189	145	149	154	231	240	242;...
                                               165	175	183	191	143	142	156	232	243	244;...
                                               167	177	182	130	136	144	158	233	245	246;...
                                               169	176	180	129	139	148	160	234	247	248;...
                                               173	170	186	131	141	155	225	235	249	250;...
                                               162	172	185	132	147	157	227	236	251	252;...
                                               164	179	190	135	138	159	226	237	253	254;...
                                               NaN	181	188	133	140	146	228	238	255	NaN];

            ChannelMap(1).ArrayLocations{4} = [193	NaN	200	NaN	213	219;...
                                               NaN	205	NaN	209	NaN	220;...
                                               195	NaN	202	NaN	215	NaN;...
                                               NaN	207	NaN	210	NaN	221;...
                                               197	NaN	204	NaN	216	NaN;...
                                               NaN	194	NaN	212	NaN	222;...
                                               199	NaN	206	NaN	217	NaN;...
                                               NaN	196	NaN	211	NaN	224;...
                                               201	NaN	208	NaN	218	NaN;...
                                               203	198	NaN	214	NaN	223];

            ChannelMap(1).ArrayOffsets = [0, 0,...
                                          sum(~isnan(ChannelMap(1).ArrayLocations{1}(:))),...
                                          sum(~isnan(ChannelMap(1).ArrayLocations{2}(:)))];


        case 'CRS02'
            ChannelMap(1).ArrayNames = {'LateralMotor', 'LateralSensory', 'MedialMotor', 'MedialSensory'};
            ChannelMap(1).ArrayRotations = [105, 60, 15, 15];
            ChannelMap(1).IsSensory = logical([0,1,0,1]);
            ChannelMap(1).IsMotor = logical([1,0,1,0]);
            ChannelMap(1).IsMedial = logical([0,0,1,1]);
            ChannelMap(1).IsLateral = logical([1,1,0,0]);
            ChannelMap(1).ArrayLocations{1} = [NaN	NaN	42	58	3	13	27	97	NaN	NaN;...
                                               NaN	34	44	57	4	19	29	98	107	NaN;...
                                               33	36	51	62	7	10	31	99	108	117;...
                                               35	38	53	60	5	12	18	100	109	119;...
                                               37	40	50	59	6	23	22	101	110	121;...
                                               39	43	46	64	9	25	24	102	111	123;...
                                               41	47	56	61	17	21	26	103	113	125;...
                                               45	49	55	63	15	14	28	104	112	127;...
                                               NaN	48	54	2	8	16	30	105	115	NaN;...
                                               NaN	NaN	52	1	11	20	32	106	NaN	NaN];

            ChannelMap(1).ArrayLocations{2} = [65	NaN	72	NaN	85	91;...
                                               NaN	77	NaN	81	NaN	92;...
                                               67	NaN	74	NaN	87	NaN;...
                                               NaN	NaN	NaN	82	NaN	94;...
                                               69	79	76	NaN	88	NaN;...
                                               NaN	66	NaN	84	NaN	93;...
                                               71	NaN	78	NaN	89	NaN;...
                                               NaN	68	NaN	83	NaN	96;...
                                               73	NaN	80	NaN	90	NaN;...
                                               75	70	NaN	86	NaN	95];

            ChannelMap(1).ArrayLocations{3} = [NaN	NaN	170	186	131	141	155	225	NaN	NaN;...
                                               NaN	162	172	185	132	147	157	226	235	NaN;...
                                               161	164	179	190	135	138	159	227	236	245;...
                                               163	166	181	188	133	140	146	228	237	247;...
                                               165	168	178	187	134	151	150	229	238	249;...
                                               167	171	174	192	137	153	152	230	239	251;...
                                               169	175	184	189	145	149	154	231	241	253;...
                                               173	177	183	191	143	142	156	232	240	255;...
                                               NaN	176	182	130	136	144	158	233	243	NaN;...
                                               NaN	NaN	180	129	139	148	160	234	NaN	NaN];
            
            ChannelMap(1).ArrayLocations{4} = [193	NaN	200	NaN	213	219;...
                                               NaN	205	NaN	209	NaN	220;...
                                               195	NaN	202	NaN	215	NaN;...
                                               NaN	207	NaN	210	NaN	222;...
                                               197	NaN	204	NaN	216	NaN;...
                                               NaN	194	NaN	212	NaN	221;...
                                               199	NaN	206	NaN	217	NaN;...
                                               NaN	196	NaN	211	NaN	224;...
                                               201	NaN	208	NaN	218	NaN;...
                                               203	198	NaN	214	NaN	223];
            ChannelMap(1).ArrayOffsets = [0, 0,...
                                          sum(~isnan(ChannelMap(1).ArrayLocations{1}(:))),...
                                          sum(~isnan(ChannelMap(1).ArrayLocations{2}(:)))];

    end
    
    % Resort the array numbers
    for i = 1:length(ChannelMap(1).ArrayLocations)
        [sorted_locs, sorted_idx] = sort(ChannelMap(1).ArrayLocations{i}(:));
        [~,resort_idx] = sort(sorted_idx);
        ChannelNumbers = [1:sum(~isnan(sorted_locs)), NaN(1,sum(isnan(sorted_locs)))] + ChannelMap(1).ArrayOffsets(i);
        ChannelMap(1).ChannelNumbers{i} = reshape(ChannelNumbers(resort_idx), size(ChannelMap(1).ArrayLocations{i}));
    end
end
