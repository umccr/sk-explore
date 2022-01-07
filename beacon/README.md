# beacon-skanwal
A federated ecosystem that allows sharing of genomic and clinical data


## Notes

The serverless beacon is up and running at [beacon](https://avt0y2mcq9.execute-api.ap-southeast-2.amazonaws.com/prod) (courtesy Florian).

### Data

Curently, it is pointing to one of the AGHA [gVCFs](s3://agha-gdr-store-dev/mitochondrial_disease/17W001163.dedup.realigned.recalibrated.hc.gvcf.gz) 

**Load data into beacon**

```
curl --header "Content-Type: application/json" \
  --request POST \
  --data '{"id":"t1","name":"test1","assemblyId":"GRCh37","vcfLocations": ["s3://umccr-primary-data-dev/bcbio_test_project_new/final/2018-12-05_2018-11-30T1125_Avner_WGS-merged/2016_249_18_WH_P017_1-ensemble-annotated.vcf.gz"]}' \
  https://avt0y2mcq9.execute-api.ap-southeast-2.amazonaws.com/prod/submit
```

**Observations**

Tried loading more data (one of the Avner sample) to the serverless beacon using:

```
curl --header "Content-Type: application/json" \
--request POST \
--data '{"id":"t2","name":"test2","assemblyId":"GRCh37","vcfLocations": ["s3://umccr-primary-data-prod/Avner/2018-12-05/final/CCR180136_WH18F001P025/2016_249_18_WH_P025_1-sv-prioritize-manta.vcf.gz"]}' \
https://avt0y2mcq9.execute-api.ap-southeast-2.amazonaws.com/prod/submit
```

The command did not return any messgae (success/error), which might be an issue as only after querying, the user can find out if the data is successfully loaded. 

And to copy data from dev bucket, used:

```
curl --header "Content-Type: application/json" \
> --request POST \
> --data '{"id":"t3","name":"test3","assemblyId":"GRCh37","vcfLocations": ["s3://umccr-primary-data-dev/test_projects/beacon_test_data/2016_249_18_WH_P025_1-sv-prioritize-manta.vcf.gz"]}' \
> https://avt0y2mcq9.execute-api.ap-southeast-2.amazonaws.com/prod/submit
```

### Example query

```
curl "https://avt0y2mcq9.execute-api.ap-southeast-2.amazonaws.com/prod/query?assemblyId=hg19;referenceName=1;referenceBases=C;alternateBases=T;start=10491;end=10491;includeDatasetResponses=HIT"
```

**Range queries**

- Deletion which starts somewhere between positions 5146261 and 5146265, and ends between 5146262 and 5374408, on chromosome 1. 

```curl "https://74palqdmj8.execute-api.ap-southeast-2.amazonaws.com/prod/query?assemblyId=GRCh37;referenceName=1;referenceBases=AT;startMin=5146261;startMax=5146265;endMin=5146262;endMax=5374408;variantType=DEL;includeDatasetResponses=HIT"```

- In the current test vcf (`2016_249_18_WH_P016_2-ensemble-annotated.vcf.gz`), there are a couple of deletions at positions `5146262` and `5374408` on chromosome 1.

To check if it's possible to pickup these two deletions in the test dataset by specifying a range value and `N` as a reference, I tested the following query:

```
curl "https://74palqdmj8.execute-api.ap-southeast-2.amazonaws.com/prod/query?assemblyId=GRCh37;referenceName=1;referenceBases=N;startMin=5146260;startMax=5374411;endMin=5146261;endMax=5374414;variantType=DEL;includeDatasetResponses=HIT"
```

It successfully returned the following  result:

```
{"beaconId": "org.umccr.beacon", "apiVersion": null, "alleleRequest": {"assemblyId": "GRCh37", "endMax": 5374414, "endMin": 5146261, "includeDatasetResponses": "HIT", "referenceBases": "N", "referenceName": "1", "startMax": 5374411, "startMin": 5146260, "variantType": "DEL", "datasetIds": null}, "exists": true, "datasetAlleleResponses": [{"datasetId": "t4", "exists": true, "frequency": 0.25, "variantCount": 2, "callCount": 2, "sampleCount": 2, "note": null, "externalUrl": null, "info": null, "error": null}]}
```

### Pointers

- [Beacon network](https://beacon-network.org/#/) - a global search engine for genetic mutations, provides a nice UI for searching variants across world's public beacons. In addition, it allows to limit the search to organizations of interest.   
- Peformance comparison between serverless abd traditional beacons
- Try loading more data into the serverless beacon
- Create some (interesting and meaningful) queries
- Is providing a referenceBases mandatory? I get an error, when not specifying 
	- It's possible to use `N` as a referenceBases value, specifically useful for range queries. 



