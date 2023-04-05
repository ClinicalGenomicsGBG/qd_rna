# SLIMS module for Cellophane

Fetch sample information from a SLIMS database. If samples are present (from a previous module or a samples file) these samples will be augmented with information from slims. If one or more sample IDs are specified they will be searched for. If no samples are present and no sample IDs are specified all samples matching the specified criteria will be fetched.

Sometimes the records holding sample information are derived from a parent record holding information linking the sample to a specific analysis. In this case the `derived_from.*` parameters can be used to initially search for parent objects that sample objects need to be derived from.

Optionally, the `bioinfo.*` parameters can be used to add/update/check "biofinformatics" objects to track the state of a secondary analysis. On pipeline completion/failure the specified field on the bioinformatics objects will be updated. to reflect this state. The state of these objects can then be checked by the wrapper to determine what samples should be excluded from analysis.

While it is possible to manually use the [mixin](#Mixins) methods manually, it is recomended to rely on the provided hooks.

## Configuration

Option                            | Type      | Required | Default | Description
----------------------------------|-----------|----------|---------|-------------
`slims.url`                       | str       | x        |         | SLIMS Server URL
`slims.username`                  | str       | x        |         | SLIMS username
`slims.password`                  | str       | x        |         | SLIMS password
`slims.content_type`              | int       | x        |         | Content type PK for sample records
`slims.map`                       | list[str] |          |         | Mapping of keys to SLIMS field(s) (see [Fields](#Fields)/[Mappings](#Mappings))
`slims.derive`                    | list[str] |          |         | Mapping of keys to SLIMS field(s) (see [Derivations](#Derivations)/[Fields](#Fields)/[Mappings](#Mappings))
`slims.find_criteria`             | str       | x        |         | SLIMS criteria for finding records (see [Criteria](#Criteria))
`slims.check_criteria`            | str       |          |         | SLIMS criteria for checking completed records (see [Criteria](#Criteria))
`slims.id`                        | list[str] |          |         | Manually select SLIMS Sample ID(s)
`slims.allow_duplicates`          | bool      |          | false   | Allow duplicate samples (eg. if a pre-hook can handle this)
`slims.dry_run`                   | bool      |          | false   | Do not create/update SLIMS bioinformatics objects
`slims.novel_max_age`             | str       |          | 1 year  | Maximum age of novel samples (eg. "4 days", "2 months", "1 year")
`slims.unrestrict_parents`        | bool      |          | false   | Allow derived records to have different IDs than parent records


### Example

```yaml
slims:
  # Ask your IT people for credentials (Remember to say please)
  url: 'http://slims.example.com/slimsrest'
  username: apiusername
  password: apipassword
  # Find sample records with content type 22 (Fastq) and qdrna in the samplesheet decription field
  # Exclude samples explicitly marked as "Do not include"
  find_criteria: |
    cntn_fk_contentType equals 22
    and cntn_cstm_doNotInclude not_equals true
    and cntn_cstm_rawSheetDescription contains qdrna
  # Check for derived records (->) with contentType 23 (Bioinformatics) and complete secondary analysis state
  check_criteria: |
    -> cntn_fk_contentType equals 23
    and cntn_cstm_SecondaryAnalysisState equals complete
  # Get file path from slims (unless specified in samples file)
  # Add backup remote keys for HCP hook
  # Add run information to samples
  map:
    - files=json:cntn_cstm_demuxerSampleResult.fastq_paths
      run=cntn_cstm_runTag
      backup=json:cntn_cstm_demuxerBackupSampleResult.remote_keys
  # Create derived records with content type = 23 (Bioinformatics)
  # Include state from cellophane sample (access key using curly braces)
  # Other parameters are required by SLIMS or records will not be created
  derive:
    - cntn_fk_contentType=23
      cntn_fk_location=83
      cntn_cstm_SecondaryAnalysisState={state}
      cntn_status=10
```


## Hooks

Name                   | When | Condition | Description
-----------------------|------|-----------|-------------
`slims_fetch`          | Pre  |           | Fetch sample info from SLIMS 
`slims_bioinformatics` | Pre  |           | Add bioinformatics records for samples
`slims_update`         | Post | Always    | Update bioinformatics state in SLIMS

## Mixins

`SlimsSamples`

```
from_records(records: list[Record], config: cfg.Config) -> data.Samples
```
Create a new data.Samples object from a list of SLIMS records.

```
update_derived(config: cfg.Config) -> None
```
Update bioinformatics state for all samples in SLIMS.

`SlimsSample`

```
from_record(record: Record, config: cfg.Config) -> data.Samples
```
Create a new data.Sample object from a single SLIMS record.

```
update_derived(str, config: cfg.Config) -> None
```
Update bioinformatics state for single sample in SLIMS.


## Criteria

Criteria are specified using `<FIELD> <OPERATOR> <VALUE>` syntax where operators take 1, 2, or more values.

> **WARNING** Invalid field names are silently ignored. This makes it very easy to accidentally fetch more samples than expected. This behavior is inherent to the official SLIMS API python bindings.

The following operators are supported:

- `equals` / `equals_ignore_case` / `not_equals` / `not_equals_ignore_case` (1 value)
- `one_of` / `not_one_of` (n values)
- `contains` / `not_contains` (n values)
- `starts_with` / `ends_with` / `not_starts_with` / `not_ends_with` (1 value)
- `between` / `not_between` (2 values)
- `greater_than` / `less_than` (1 value)

Complex boolean criteria can be constructed using `and`/`or` and parentheses.

To find derived records, the "->" operator can be used between two criteria. This will find all records that match the first criteria and then find all records that match the second criteria and are derived from the first (and so on). This operator can not be used inside parentheses. In order to limit the search space, parent records are expected to have the same ID as the derived records. This behaviour can be disabled using the `slims.unrestrict_parents` option. The `slims.bioinfo_check_criteria` currently does not support this operator.

eg. `cntn_cstm_foo equals a -> cntn_cstm_foo equals b and (cntn_cstm_bar not_between_inclusive c d or cntn_cstm_baz equals e)`

## Fields

Fields are specified using the full SLIMS field name (eg. `cntn_fk_contentType`). Sub-fields inside JSON fields can be accessed by using the `json:` prefix and dot-notation (eg. json:cntn_cstm_foo.bar.baz) 

## Mappings

SLIMS fields can be mapped to keys on the `data.Samples` object using `<KEY>=<FIELD>` with the same field syntax as above. This can be used to override any field on the `data.Samples` object (eg. `files=json:cntn_cstm_myJSONField.fastq_paths`). Use whitespace/newline to separate each `<KEY>=<FIELD>` pair in the desired mapping.

## Derivations

Derived records can be automatically created by specifying a mapping using `<FIELD>=<VALUE>` syntax (eg. `cntn_fk_contentType=1234 cntn_cstm_foo={state}`). Curly braces can be used to reference keys on the Sample object. The value of the field will be used as the value for the derived record. The derived record will be created if it does not exist and updated if it does. The derived record will be linked to the parent record using the `cntn_fk_originalContent` field. The `cntn_fk_contentType` field must be specified on the derived record. SLIMS may fail to create a records silently if required fields are not specified, and there is currently no way to explicitly check for required fields.

To create multiple derived records, multiple mappings can be specified as a list (or by passing the `--slims_derive` flag multiple times).