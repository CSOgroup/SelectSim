# Filter a MAF dataframe by a combination of column values

Filter a MAF dataframe by a combination of column values

## Usage

``` r
filter_maf_complex(maf, values, ...)
```

## Arguments

- maf:

  A MAF dataframe

- values:

  A dataframe of (column, value) pairs to match against

- ...:

  Additional arguments passed to merge

## Value

Filtered MAF dataframe containing only rows matching the value
combinations.

## Examples

``` r
data(luad_maf, package = "SelectSim")
combos <- data.frame(Hugo_Symbol = "TP53",
                     Variant_Classification = "Missense_Mutation")
filter_maf_complex(luad_maf, combos,
                   by.x = c("Hugo_Symbol","Variant_Classification"),
                   by.y = c("Hugo_Symbol","Variant_Classification"))
#>     Hugo_Symbol Variant_Classification Chromosome Start_Position End_Position
#> 1          TP53      Missense_Mutation         17        7577547      7577547
#> 2          TP53      Missense_Mutation         17        7577547      7577547
#> 3          TP53      Missense_Mutation         17        7578507      7578507
#> 4          TP53      Missense_Mutation         17        7579373      7579373
#> 5          TP53      Missense_Mutation         17        7577574      7577574
#> 6          TP53      Missense_Mutation         17        7577108      7577108
#> 7          TP53      Missense_Mutation         17        7577121      7577121
#> 8          TP53      Missense_Mutation         17        7578475      7578475
#> 9          TP53      Missense_Mutation         17        7577097      7577097
#> 10         TP53      Missense_Mutation         17        7577507      7577507
#> 11         TP53      Missense_Mutation         17        7574017      7574017
#> 12         TP53      Missense_Mutation         17        7579361      7579361
#> 13         TP53      Missense_Mutation         17        7577117      7577117
#> 14         TP53      Missense_Mutation         17        7578461      7578461
#> 15         TP53      Missense_Mutation         17        7578403      7578403
#> 16         TP53      Missense_Mutation         17        7577551      7577551
#> 17         TP53      Missense_Mutation         17        7578419      7578419
#> 18         TP53      Missense_Mutation         17        7578455      7578455
#> 19         TP53      Missense_Mutation         17        7577551      7577551
#> 20         TP53      Missense_Mutation         17        7578475      7578475
#> 21         TP53      Missense_Mutation         17        7577100      7577100
#> 22         TP53      Missense_Mutation         17        7578457      7578457
#> 23         TP53      Missense_Mutation         17        7578266      7578266
#> 24         TP53      Missense_Mutation         17        7578550      7578550
#> 25         TP53      Missense_Mutation         17        7578542      7578542
#> 26         TP53      Missense_Mutation         17        7578536      7578536
#> 27         TP53      Missense_Mutation         17        7577096      7577096
#> 28         TP53      Missense_Mutation         17        7577547      7577547
#> 29         TP53      Missense_Mutation         17        7578448      7578448
#> 30         TP53      Missense_Mutation         17        7577509      7577509
#> 31         TP53      Missense_Mutation         17        7577124      7577124
#> 32         TP53      Missense_Mutation         17        7577100      7577100
#> 33         TP53      Missense_Mutation         17        7577082      7577082
#> 34         TP53      Missense_Mutation         17        7578466      7578466
#> 35         TP53      Missense_Mutation         17        7578236      7578236
#> 36         TP53      Missense_Mutation         17        7577082      7577082
#> 37         TP53      Missense_Mutation         17        7577538      7577538
#> 38         TP53      Missense_Mutation         17        7578271      7578271
#> 39         TP53      Missense_Mutation         17        7578427      7578427
#> 41         TP53      Missense_Mutation         17        7579358      7579358
#> 42         TP53      Missense_Mutation         17        7578268      7578268
#> 43         TP53      Missense_Mutation         17        7578413      7578413
#> 44         TP53      Missense_Mutation         17        7578454      7578454
#> 45         TP53      Missense_Mutation         17        7574026      7574026
#> 46         TP53      Missense_Mutation         17        7577120      7577120
#> 47         TP53      Missense_Mutation         17        7577114      7577114
#> 48         TP53      Missense_Mutation         17        7577548      7577548
#> 49         TP53      Missense_Mutation         17        7578190      7578190
#> 50         TP53      Missense_Mutation         17        7577121      7577121
#> 51         TP53      Missense_Mutation         17        7578406      7578406
#> 52         TP53      Missense_Mutation         17        7578443      7578443
#> 53         TP53      Missense_Mutation         17        7577566      7577566
#> 54         TP53      Missense_Mutation         17        7578403      7578403
#> 55         TP53      Missense_Mutation         17        7577095      7577095
#> 56         TP53      Missense_Mutation         17        7577535      7577535
#> 58         TP53      Missense_Mutation         17        7577141      7577141
#> 59         TP53      Missense_Mutation         17        7578461      7578461
#> 60         TP53      Missense_Mutation         17        7578403      7578403
#> 61         TP53      Missense_Mutation         17        7578526      7578526
#> 62         TP53      Missense_Mutation         17        7578455      7578455
#> 63         TP53      Missense_Mutation         17        7577574      7577574
#> 64         TP53      Missense_Mutation         17        7577556      7577556
#> 65         TP53      Missense_Mutation         17        7578442      7578442
#> 66         TP53      Missense_Mutation         17        7578203      7578203
#> 67         TP53      Missense_Mutation         17        7578271      7578271
#> 68         TP53      Missense_Mutation         17        7578406      7578406
#> 69         TP53      Missense_Mutation         17        7577123      7577123
#> 70         TP53      Missense_Mutation         17        7577117      7577117
#> 71         TP53      Missense_Mutation         17        7578498      7578498
#> 72         TP53      Missense_Mutation         17        7577587      7577587
#> 73         TP53      Missense_Mutation         17        7577085      7577085
#> 74         TP53      Missense_Mutation         17        7576853      7576853
#> 75         TP53      Missense_Mutation         17        7577556      7577556
#> 76         TP53      Missense_Mutation         17        7574027      7574027
#> 77         TP53      Missense_Mutation         17        7577120      7577120
#> 78         TP53      Missense_Mutation         17        7574017      7574017
#> 79         TP53      Missense_Mutation         17        7578455      7578455
#> 80         TP53      Missense_Mutation         17        7577565      7577565
#> 81         TP53      Missense_Mutation         17        7579374      7579374
#> 82         TP53      Missense_Mutation         17        7578457      7578457
#> 83         TP53      Missense_Mutation         17        7578448      7578448
#> 84         TP53      Missense_Mutation         17        7579373      7579373
#> 85         TP53      Missense_Mutation         17        7577536      7577536
#> 86         TP53      Missense_Mutation         17        7577535      7577535
#> 87         TP53      Missense_Mutation         17        7578235      7578235
#> 88         TP53      Missense_Mutation         17        7578534      7578534
#> 91         TP53      Missense_Mutation         17        7578457      7578457
#> 92         TP53      Missense_Mutation         17        7577120      7577120
#> 93         TP53      Missense_Mutation         17        7578443      7578443
#> 94         TP53      Missense_Mutation         17        7577570      7577570
#> 95         TP53      Missense_Mutation         17        7578401      7578401
#> 96         TP53      Missense_Mutation         17        7577099      7577099
#> 97         TP53      Missense_Mutation         17        7574018      7574018
#> 98         TP53      Missense_Mutation         17        7577541      7577541
#> 99         TP53      Missense_Mutation         17        7577120      7577120
#> 100        TP53      Missense_Mutation         17        7577536      7577536
#> 101        TP53      Missense_Mutation         17        7578461      7578461
#> 103        TP53      Missense_Mutation         17        7577535      7577535
#> 104        TP53      Missense_Mutation         17        7578550      7578550
#> 105        TP53      Missense_Mutation         17        7578508      7578508
#> 106        TP53      Missense_Mutation         17        7578268      7578268
#> 107        TP53      Missense_Mutation         17        7577121      7577121
#> 108        TP53      Missense_Mutation         17        7577548      7577548
#> 109        TP53      Missense_Mutation         17        7577570      7577570
#> 110        TP53      Missense_Mutation         17        7578272      7578272
#> 111        TP53      Missense_Mutation         17        7577548      7577548
#> 112        TP53      Missense_Mutation         17        7573964      7573964
#> 113        TP53      Missense_Mutation         17        7577094      7577094
#> 114        TP53      Missense_Mutation         17        7577081      7577081
#> 115        TP53      Missense_Mutation         17        7577547      7577547
#> 116        TP53      Missense_Mutation         17        7574018      7574018
#> 117        TP53      Missense_Mutation         17        7578442      7578442
#> 118        TP53      Missense_Mutation         17        7578457      7578457
#> 119        TP53      Missense_Mutation         17        7574026      7574026
#> 120        TP53      Missense_Mutation         17        7577538      7577538
#> 121        TP53      Missense_Mutation         17        7577082      7577082
#> 122        TP53      Missense_Mutation         17        7577127      7577127
#> 123        TP53      Missense_Mutation         17        7577534      7577534
#> 124        TP53      Missense_Mutation         17        7578457      7578457
#> 125        TP53      Missense_Mutation         17        7578395      7578395
#> 126        TP53      Missense_Mutation         17        7578208      7578208
#> 127        TP53      Missense_Mutation         17        7577099      7577099
#> 128        TP53      Missense_Mutation         17        7577538      7577538
#> 129        TP53      Missense_Mutation         17        7578235      7578235
#> 130        TP53      Missense_Mutation         17        7578264      7578264
#> 131        TP53      Missense_Mutation         17        7578461      7578461
#> 132        TP53      Missense_Mutation         17        7578392      7578392
#> 133        TP53      Missense_Mutation         17        7578271      7578271
#> 134        TP53      Missense_Mutation         17        7579518      7579518
#> 135        TP53      Missense_Mutation         17        7577097      7577097
#> 136        TP53      Missense_Mutation         17        7577085      7577085
#> 137        TP53      Missense_Mutation         17        7574017      7574017
#> 138        TP53      Missense_Mutation         17        7577539      7577539
#> 139        TP53      Missense_Mutation         17        7577559      7577559
#> 140        TP53      Missense_Mutation         17        7577547      7577547
#> 141        TP53      Missense_Mutation         17        7577120      7577120
#> 143        TP53      Missense_Mutation         17        7579358      7579358
#> 144        TP53      Missense_Mutation         17        7578499      7578499
#> 146        TP53      Missense_Mutation         17        7578393      7578393
#> 147        TP53      Missense_Mutation         17        7578190      7578190
#> 148        TP53      Missense_Mutation         17        7578457      7578457
#> 149        TP53      Missense_Mutation         17        7578457      7578457
#> 150        TP53      Missense_Mutation         17        7577120      7577120
#> 151        TP53      Missense_Mutation         17        7577090      7577090
#> 152        TP53      Missense_Mutation         17        7579329      7579329
#> 153        TP53      Missense_Mutation         17        7578271      7578271
#> 154        TP53      Missense_Mutation         17        7578235      7578235
#> 155        TP53      Missense_Mutation         17        7577153      7577153
#> 156        TP53      Missense_Mutation         17        7577130      7577130
#> 157        TP53      Missense_Mutation         17        7577099      7577099
#> 158        TP53      Missense_Mutation         17        7577141      7577141
#> 159        TP53      Missense_Mutation         17        7577097      7577097
#> 160        TP53      Missense_Mutation         17        7578457      7578457
#> 162        TP53      Missense_Mutation         17        7577082      7577082
#> 163        TP53      Missense_Mutation         17        7579457      7579457
#> 164        TP53      Missense_Mutation         17        7577108      7577108
#> 165        TP53      Missense_Mutation         17        7578406      7578406
#> 166        TP53      Missense_Mutation         17        7577114      7577114
#> 167        TP53      Missense_Mutation         17        7578469      7578469
#> 168        TP53      Missense_Mutation         17        7578235      7578235
#> 169        TP53      Missense_Mutation         17        7577535      7577535
#> 171        TP53      Missense_Mutation         17        7577538      7577538
#> 172        TP53      Missense_Mutation         17        7577532      7577532
#> 173        TP53      Missense_Mutation         17        7577556      7577556
#> 174        TP53      Missense_Mutation         17        7577094      7577094
#> 176        TP53      Missense_Mutation         17        7577129      7577129
#> 177        TP53      Missense_Mutation         17        7578463      7578463
#> 178        TP53      Missense_Mutation         17        7577105      7577105
#> 179        TP53      Missense_Mutation         17        7577568      7577568
#> 180        TP53      Missense_Mutation         17        7577138      7577138
#>             Tumor_Sample_Barcode          sample HGVSp_Short
#> 1   TCGA-97-7937-01A-11D-2167-08 TCGA-97-7937-01     p.G245V
#> 2   TCGA-55-8301-01A-11D-2284-08 TCGA-55-8301-01     p.G245V
#> 3   TCGA-93-A4JO-01A-21D-A24P-08 TCGA-93-A4JO-01     p.C141W
#> 4   TCGA-95-8039-01A-11D-2238-08 TCGA-95-8039-01     p.G105D
#> 5   TCGA-MP-A4SV-01A-11D-A24P-08 TCGA-MP-A4SV-01     p.Y236C
#> 6   TCGA-49-4505-01A-01D-1931-08 TCGA-49-4505-01     p.C277F
#> 7   TCGA-78-7540-01A-11D-2063-08 TCGA-78-7540-01     p.R273C
#> 8   TCGA-L4-A4E5-01A-11D-A24P-08 TCGA-L4-A4E5-01     p.P152L
#> 9   TCGA-69-A59K-01A-11D-A25L-08 TCGA-69-A59K-01     p.D281Y
#> 10  TCGA-86-6851-01A-11D-1945-08 TCGA-86-6851-01     p.E258D
#> 11  TCGA-44-7670-01A-11D-2063-08 TCGA-44-7670-01     p.R337L
#> 12  TCGA-44-5644-01A-21D-2036-08 TCGA-44-5644-01     p.F109C
#> 13  TCGA-62-8394-01A-11D-2323-08 TCGA-62-8394-01     p.V274G
#> 14  TCGA-53-A4EZ-01A-12D-A24P-08 TCGA-53-A4EZ-01     p.V157F
#> 15  TCGA-73-4676-01A-01D-1753-08 TCGA-73-4676-01     p.C176F
#> 16  TCGA-78-7220-01A-11D-2036-08 TCGA-78-7220-01     p.G244C
#> 17  TCGA-86-8056-01A-11D-2238-08 TCGA-86-8056-01     p.E171K
#> 18  TCGA-78-7145-01A-11D-2036-08 TCGA-78-7145-01     p.A159P
#> 19  TCGA-55-1592-01A-01D-0969-08 TCGA-55-1592-01     p.G244C
#> 20  TCGA-L9-A444-01A-21D-A24D-08 TCGA-L9-A444-01     p.P152L
#> 21  TCGA-64-5775-01A-01D-1625-08 TCGA-64-5775-01     p.R280G
#> 22  TCGA-49-AARQ-01A-11D-A410-08 TCGA-49-AARQ-01     p.R158L
#> 23  TCGA-L9-A7SV-01A-11D-A397-08 TCGA-L9-A7SV-01     p.I195F
#> 24  TCGA-44-A479-01A-31D-A24D-08 TCGA-44-A479-01     p.S127F
#> 25  TCGA-75-5125-01A-01D-1753-08 TCGA-75-5125-01     p.L130F
#> 26  TCGA-44-6779-01A-11D-1855-08 TCGA-44-6779-01     p.K132E
#> 27  TCGA-05-4410-01A-21D-1855-08 TCGA-05-4410-01     p.D281V
#> 28  TCGA-44-6777-01A-11D-1855-08 TCGA-44-6777-01     p.G245V
#> 29  TCGA-MP-A4TI-01A-21D-A24P-08 TCGA-MP-A4TI-01     p.A161D
#> 30  TCGA-L9-A743-01A-43D-A397-08 TCGA-L9-A743-01     p.E258Q
#> 31  TCGA-50-6594-01A-11D-1753-08 TCGA-50-6594-01     p.V272M
#> 32  TCGA-55-8620-01A-11D-2393-08 TCGA-55-8620-01     p.R280G
#> 33  TCGA-86-A4P7-01A-11D-A24P-08 TCGA-86-A4P7-01     p.E286K
#> 34  TCGA-05-4426-01A-01D-1265-08 TCGA-05-4426-01     p.T155I
#> 35  TCGA-95-7567-01A-11D-2063-08 TCGA-95-7567-01     p.Y205H
#> 36  TCGA-MN-A4N1-01A-11D-A24P-08 TCGA-MN-A4N1-01     p.E286Q
#> 37  TCGA-86-8673-01A-11D-2393-08 TCGA-86-8673-01     p.R248L
#> 38  TCGA-86-8280-01A-11D-2284-08 TCGA-86-8280-01     p.H193P
#> 39  TCGA-44-6775-01A-31D-A27T-08 TCGA-44-6775-01     p.H168L
#> 41  TCGA-91-6847-01A-11D-1945-08 TCGA-91-6847-01     p.R110L
#> 42  TCGA-62-A46Y-01A-11D-A24D-08 TCGA-62-A46Y-01     p.L194R
#> 43  TCGA-99-8025-01A-11D-2238-08 TCGA-99-8025-01     p.V173L
#> 44  TCGA-MP-A4TC-01A-11D-A24P-08 TCGA-MP-A4TC-01     p.A159V
#> 45  TCGA-55-6985-01A-11D-1945-08 TCGA-55-6985-01     p.G334V
#> 46  TCGA-55-6969-01A-11D-1945-08 TCGA-55-6969-01     p.R273L
#> 47  TCGA-83-5908-01A-21D-2284-08 TCGA-83-5908-01     p.C275F
#> 48  TCGA-44-A4SS-01A-11D-A24P-08 TCGA-44-A4SS-01     p.G245C
#> 49  TCGA-86-8074-01A-11D-2238-08 TCGA-86-8074-01     p.Y220C
#> 50  TCGA-55-A493-01A-11D-A24D-08 TCGA-55-A493-01     p.R273S
#> 51  TCGA-50-6595-01A-12D-1855-08 TCGA-50-6595-01     p.R175H
#> 52  TCGA-95-7944-01A-11D-2184-08 TCGA-95-7944-01     p.Y163D
#> 53  TCGA-J2-A4AD-01A-11D-A24D-08 TCGA-J2-A4AD-01     p.N239D
#> 54  TCGA-05-4395-01A-01D-1265-08 TCGA-05-4395-01     p.C176Y
#> 55  TCGA-05-4432-01A-01D-1265-08 TCGA-05-4432-01     p.D281E
#> 56  TCGA-35-4122-01A-01D-1105-08 TCGA-35-4122-01     p.R249M
#> 58  TCGA-55-8614-01A-11D-2393-08 TCGA-55-8614-01     p.G266V
#> 59  TCGA-78-8640-01A-11D-2393-08 TCGA-78-8640-01     p.V157L
#> 60  TCGA-86-7711-01A-11D-2063-08 TCGA-86-7711-01     p.C176F
#> 61  TCGA-55-7573-01A-11D-2036-08 TCGA-55-7573-01     p.C135F
#> 62  TCGA-38-4632-01A-01D-1753-08 TCGA-38-4632-01     p.A159P
#> 63  TCGA-62-A471-01A-12D-A24D-08 TCGA-62-A471-01     p.Y236C
#> 64  TCGA-99-8028-01A-11D-2238-08 TCGA-99-8028-01     p.C242F
#> 65  TCGA-55-7570-01A-11D-2036-08 TCGA-55-7570-01     p.Y163C
#> 66  TCGA-50-5068-01A-01D-1625-08 TCGA-50-5068-01     p.V216M
#> 67  TCGA-55-A491-01A-11D-A24D-08 TCGA-55-A491-01     p.H193R
#> 68  TCGA-49-4514-01A-21D-1855-08 TCGA-49-4514-01     p.R175H
#> 69  TCGA-05-4382-01A-01D-1931-08 TCGA-05-4382-01     p.V272G
#> 70  TCGA-MN-A4N4-01A-12D-A24P-08 TCGA-MN-A4N4-01     p.V274D
#> 71  TCGA-55-7911-01A-11D-2167-08 TCGA-55-7911-01     p.Q144H
#> 72  TCGA-44-6778-01A-11D-1855-08 TCGA-44-6778-01     p.I232F
#> 73  TCGA-50-5939-01A-11D-1625-08 TCGA-50-5939-01     p.E285K
#> 74  TCGA-62-8399-01A-21D-2323-08 TCGA-62-8399-01     p.Q331H
#> 75  TCGA-MP-A4TK-01A-11D-A24P-08 TCGA-MP-A4TK-01     p.C242F
#> 76  TCGA-69-7974-01A-11D-2184-08 TCGA-69-7974-01     p.G334W
#> 77  TCGA-55-7994-01A-11D-2184-08 TCGA-55-7994-01     p.R273L
#> 78  TCGA-55-8514-01A-11D-2393-08 TCGA-55-8514-01     p.R337P
#> 79  TCGA-64-1677-01A-01W-0928-08 TCGA-64-1677-01     p.A159P
#> 80  TCGA-64-5778-01A-01D-1625-08 TCGA-64-5778-01     p.N239S
#> 81  TCGA-49-AAR4-01A-12D-A410-08 TCGA-49-AAR4-01     p.G105C
#> 82  TCGA-95-7043-01A-11D-1945-08 TCGA-95-7043-01     p.R158L
#> 83  TCGA-69-7980-01A-11D-2184-08 TCGA-69-7980-01     p.A161V
#> 84  TCGA-05-4398-01A-01D-1265-08 TCGA-05-4398-01     p.G105D
#> 85  TCGA-50-6592-01A-11D-1753-08 TCGA-50-6592-01     p.R249G
#> 86  TCGA-05-5423-01A-01D-1625-08 TCGA-05-5423-01     p.R249M
#> 87  TCGA-55-8205-01A-11D-2238-08 TCGA-55-8205-01     p.Y205C
#> 88  TCGA-95-8494-01A-11D-2323-08 TCGA-95-8494-01     p.K132N
#> 91  TCGA-55-7907-01A-11D-2167-08 TCGA-55-7907-01     p.R158P
#> 92  TCGA-97-8179-01A-11D-2284-08 TCGA-97-8179-01     p.R273L
#> 93  TCGA-86-8055-01A-11D-2238-08 TCGA-86-8055-01     p.Y163H
#> 94  TCGA-53-7626-01A-12D-2063-08 TCGA-53-7626-01     p.M237I
#> 95  TCGA-62-8402-01A-11D-2323-08 TCGA-62-8402-01     p.P177T
#> 96  TCGA-05-4420-01A-01D-1265-08 TCGA-05-4420-01     p.R280I
#> 97  TCGA-50-6673-01A-11D-1945-08 TCGA-50-6673-01     p.R337C
#> 98  TCGA-86-8358-01A-11D-2323-08 TCGA-86-8358-01     p.N247I
#> 99  TCGA-L9-A7SV-01A-11D-A397-08 TCGA-L9-A7SV-01     p.R273L
#> 100 TCGA-97-7554-01A-11D-2036-08 TCGA-97-7554-01     p.R249W
#> 101 TCGA-62-A46U-01A-11D-A24D-08 TCGA-62-A46U-01     p.V157L
#> 103 TCGA-05-4430-01A-02D-1265-08 TCGA-05-4430-01     p.R249M
#> 104 TCGA-67-3771-01A-01D-1040-01 TCGA-67-3771-01     p.S127C
#> 105 TCGA-55-6979-01A-11D-1945-08 TCGA-55-6979-01     p.C141F
#> 106 TCGA-05-4424-01A-22D-1855-08 TCGA-05-4424-01     p.L194R
#> 107 TCGA-49-AARE-01A-11D-A410-08 TCGA-49-AARE-01     p.R273G
#> 108 TCGA-55-8301-01A-11D-2284-08 TCGA-55-8301-01     p.G245S
#> 109 TCGA-49-AARN-01A-21D-A410-08 TCGA-49-AARN-01     p.M237I
#> 110 TCGA-97-8175-01A-11D-2284-08 TCGA-97-8175-01     p.H193N
#> 111 TCGA-44-6774-01A-21D-1855-08 TCGA-44-6774-01     p.G245C
#> 112 TCGA-55-7903-01A-11D-2167-08 TCGA-55-7903-01     p.A355T
#> 113 TCGA-50-8460-01A-11D-2323-08 TCGA-50-8460-01     p.R282W
#> 114 TCGA-44-3396-01A-01D-1553-08 TCGA-44-3396-01     p.E286G
#> 115 TCGA-44-A4SS-01A-11D-A24P-08 TCGA-44-A4SS-01     p.G245V
#> 116 TCGA-55-A48Z-01A-12D-A24P-08 TCGA-55-A48Z-01     p.R337C
#> 117 TCGA-75-5147-01A-01D-1625-08 TCGA-75-5147-01     p.Y163C
#> 118 TCGA-55-8092-01A-11D-2238-08 TCGA-55-8092-01     p.R158L
#> 119 TCGA-91-6831-01A-11D-1855-08 TCGA-91-6831-01     p.G334V
#> 120 TCGA-78-7149-01A-11D-2036-08 TCGA-78-7149-01     p.R248L
#> 121 TCGA-93-A4JN-01A-11D-A24P-08 TCGA-93-A4JN-01     p.E286Q
#> 122 TCGA-78-7154-01A-11D-2036-08 TCGA-78-7154-01     p.E271K
#> 123 TCGA-55-8089-01A-11D-2238-08 TCGA-55-8089-01     p.R249S
#> 124 TCGA-50-5933-01A-11D-1753-08 TCGA-50-5933-01     p.R158L
#> 125 TCGA-55-7995-01A-11D-2184-08 TCGA-55-7995-01     p.H179Y
#> 126 TCGA-4B-A93V-01A-11D-A397-08 TCGA-4B-A93V-01     p.H214R
#> 127 TCGA-44-2662-01A-01D-A271-08 TCGA-44-2662-01     p.R280T
#> 128 TCGA-73-4666-01A-01D-1265-08 TCGA-73-4666-01     p.R248P
#> 129 TCGA-95-7947-01A-11D-2184-08 TCGA-95-7947-01     p.Y205C
#> 130 TCGA-05-5423-01A-01D-1625-08 TCGA-05-5423-01     p.I195M
#> 131 TCGA-75-6214-01A-41D-1945-08 TCGA-75-6214-01     p.V157F
#> 132 TCGA-64-1681-01A-11D-2063-08 TCGA-64-1681-01     p.E180K
#> 133 TCGA-44-A47G-01A-21D-A24D-08 TCGA-44-A47G-01     p.H193L
#> 134 TCGA-J2-A4AD-01A-11D-A24D-08 TCGA-J2-A4AD-01      p.D57N
#> 135 TCGA-55-1596-01A-01D-1040-01 TCGA-55-1596-01     p.D281Y
#> 136 TCGA-93-8067-01A-11D-2284-08 TCGA-93-8067-01     p.E285K
#> 137 TCGA-78-7155-01A-11D-2036-08 TCGA-78-7155-01     p.R337L
#> 138 TCGA-L9-A8F4-01A-11D-A397-08 TCGA-L9-A8F4-01     p.R248W
#> 139 TCGA-35-5375-01A-01D-1625-08 TCGA-35-5375-01     p.S241F
#> 140 TCGA-64-1679-01A-21D-2063-08 TCGA-64-1679-01     p.G245V
#> 141 TCGA-49-AAQV-01A-11D-A397-08 TCGA-49-AAQV-01     p.R273H
#> 143 TCGA-55-A48X-01A-11D-A24D-08 TCGA-55-A48X-01     p.R110L
#> 144 TCGA-55-6979-01A-11D-1945-08 TCGA-55-6979-01     p.Q144L
#> 146 TCGA-50-5066-01A-01D-1625-08 TCGA-50-5066-01     p.H179Q
#> 147 TCGA-91-6840-01A-11D-1945-08 TCGA-91-6840-01     p.Y220C
#> 148 TCGA-55-A4DG-01A-11D-A24D-08 TCGA-55-A4DG-01     p.R158L
#> 149 TCGA-MP-A4T4-01A-11D-A25L-08 TCGA-MP-A4T4-01     p.R158P
#> 150 TCGA-95-7039-01A-11D-1945-08 TCGA-95-7039-01     p.R273L
#> 151 TCGA-38-6178-01A-11D-1753-08 TCGA-38-6178-01     p.R283P
#> 152 TCGA-55-7910-01A-11D-2167-08 TCGA-55-7910-01     p.K120E
#> 153 TCGA-86-8279-01A-11D-2284-08 TCGA-86-8279-01     p.H193R
#> 154 TCGA-49-AAR9-01A-21D-A410-08 TCGA-49-AAR9-01     p.Y205C
#> 155 TCGA-50-5930-01A-11D-1753-08 TCGA-50-5930-01     p.G262V
#> 156 TCGA-69-7760-01A-11D-2167-08 TCGA-69-7760-01     p.F270V
#> 157 TCGA-MP-A4TA-01A-21D-A24P-08 TCGA-MP-A4TA-01     p.R280I
#> 158 TCGA-55-A57B-01A-12D-A397-08 TCGA-55-A57B-01     p.G266V
#> 159 TCGA-05-4427-01A-21D-1855-08 TCGA-05-4427-01     p.D281N
#> 160 TCGA-78-7150-01A-21D-2036-08 TCGA-78-7150-01     p.R158L
#> 162 TCGA-L9-A5IP-01A-21D-A397-08 TCGA-L9-A5IP-01     p.E286K
#> 163 TCGA-49-4487-01A-21D-1855-08 TCGA-49-4487-01      p.P77L
#> 164 TCGA-05-4405-01A-21D-1855-08 TCGA-05-4405-01     p.C277F
#> 165 TCGA-05-4405-01A-21D-1855-08 TCGA-05-4405-01     p.R175H
#> 166 TCGA-55-8208-01A-11D-2238-08 TCGA-55-8208-01     p.C275F
#> 167 TCGA-55-8506-01A-11D-2393-08 TCGA-55-8506-01     p.G154V
#> 168 TCGA-05-4384-01A-01D-1753-08 TCGA-05-4384-01     p.Y205C
#> 169 TCGA-MP-A4TF-01A-11D-A25L-08 TCGA-MP-A4TF-01     p.R249M
#> 171 TCGA-86-7954-01A-11D-2184-08 TCGA-86-7954-01     p.R248L
#> 172 TCGA-78-7146-01A-11D-2036-08 TCGA-78-7146-01     p.P250L
#> 173 TCGA-97-A4LX-01A-11D-A24P-08 TCGA-97-A4LX-01     p.C242F
#> 174 TCGA-97-7553-01A-21D-2036-08 TCGA-97-7553-01     p.R282W
#> 176 TCGA-78-7536-01A-11D-2063-08 TCGA-78-7536-01     p.F270C
#> 177 TCGA-91-6848-01A-11D-1945-08 TCGA-91-6848-01     p.R156P
#> 178 TCGA-75-5146-01A-01D-1625-08 TCGA-75-5146-01     p.P278H
#> 179 TCGA-05-5428-01A-01D-1625-08 TCGA-05-5428-01     p.C238Y
#> 180 TCGA-75-5146-01A-01D-1625-08 TCGA-75-5146-01     p.R267L
```
