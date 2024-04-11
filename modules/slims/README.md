# SLIMS module for Cellophane

Fetch sample information from a SLIMS database. If samples are present (from a previous module or a samples file) these samples will be augmented with information from slims. If one or more sample IDs are specified they will be searched for. If no samples are present and no sample IDs are specified all samples matching the specified criteria will be fetched.

## Configuration

Option                            | Type      | Required | Default | Description
----------------------------------|-----------|----------|---------|-------------
`slims.url`                       | string    | Yes      |         | SLIMS URL
`slims.username`                  | string    | Yes      |         | SLIMS username
`slims.password`                  | string    | Yes      |         | SLIMS password
`slims.map`                       | mapping   |          |         | Mapping of keys to SLIMS field(s) (Use json: prefix and dot notation for JSON fields)
`slims.derive`                    | array     |          |         | Mapping for creating derived records in SLIMS (Use curly braces to access
`slims.find_criteria`             | string    |          |         | SLIMS criteria for finding records (eg. "cntn_cstm_SecondaryAnalysis equals 1337")
`slims.check_criteria`            | string    |          |         | SLIMS criteria for checking completed records (eg. "cntn_cstm_SecondaryAnalysis equals 1337")
`slims.id`                        | array     |          |         | Manually select SLIMS Sample ID(s)
`slims.allow_duplicates`          | bool      |          | False   | Allow duplicate samples (eg. if a pre-hook can handle this)
`slims.unrestrict_parents`        | bool      |          | False   | Allow parent records to have different IDs than the child records
`slims.dry_run`                   | bool      |          | False   | Do not create SLIMS bioinformatics objects
`slims.novel_max_age`             | string    |          | 1 year  | Maximum age of novel samples (eg. "7 days", "1 year")


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
    files: json:cntn_cstm_demuxerSampleResult.fastq_paths
    run: cntn_cstm_runTag
    backup: json:cntn_cstm_demuxerBackupSampleResult.remote_keys
  # Create derived records with content type = 23 (Bioinformatics)
  # Include state from cellophane sample (access key using curly braces)
  # Other parameters are required by SLIMS or records will not be created
  derive:
    cntn_fk_contentType: 23
    cntn_fk_location: 83
    cntn_cstm_SecondaryAnalysisState: '{sample.state}'
    cntn_status: 10
```

## Hooks

Name                   | When | Condition | Description
-----------------------|------|-----------|-------------
`slims_fetch`          | Pre  |           | Fetch sample info
`slims_derive`         | Pre  |           | Create derived records
`slims_running`        | Pre  |           | Mark samples as running and update derived records
`slims_update`         | Post | Always    | Update derived records

## Mixins

### Samples

```python
Samples.from_records(
  records: list[Record],
  config: cfg.Config,
) -> data.Samples
```

Classmethod to create a new data.Samples object from a list of SLIMS records.

```python
Samples.from_criteria(
  records: list[Record],
  config: cfg.Config,
  connection: slims.SLIMS | None = None,
) -> data.Samples
```

Classmethod to create a new data.Samples object with records matching the specified [criteria](criteria). Optionally pass a SLIMS connection object to avoid creating a new one.

```python
Samples.update_derived(
  config: cfg.Config,
) -> None
```

Update bioinformatics state for all samples in SLIMS.

```python
Samples.set_state(
  value: Literal["novel", "running", "complete", "failed"],
) -> None
```

Update state for all samples.

---

### Sample

```python
Sample.derived: list[Record] | None
```

Derived records for this sample.

```python
Sample.record: slims.Record | None
```

Record object for this sample (or None if sample has no Record).

```python
Sample.pk: str | None
```

SLIMS primary key for this sample (or None if sample has no Record).

```python
Sample.connection: slims.Record | None
```

Cached SLIMS connection object for this sample (or None if sample has no Record).

```python
Sample.state: Literal["novel", "running", "complete", "failed"]
```

State of the sample. Used to track sample state in SLIMS

```python
Sample.from_record(
  record: slims.Record,
  config: cellophane.Config,
) -> data.Samples
```

Create a new data.Sample object from a single SLIMS record.

```python
Sample.update_derived(
  config: cellophane.Config,
) -> None
```

Update bioinformatics state for single sample in SLIMS.

## Criteria

Criteria are specified using `<FIELD> <OPERATOR> <VALUE> ...` syntax where operators may take one or more values.

> **WARNING** Invalid field names are silently ignored. This makes it very easy to accidentally fetch more samples than expected. This behavior is inherent to the official SLIMS API python bindings.

The following operators are supported:

Operator             | Inverse                  | Value(s) | Description
---------------------|--------------------------|----------|-------------
`equals`             | `not_equals`             | 1        | Field equals (or does not equal) value
`equals_ignore_case` | `not_equals_ignore_case` | 1        | Field equals (or does not equal) value (case insensitive)
`one_of`             | `not_one_of`             | n        | Field is one of (or not one of) values
`contains`           | `not_contains`           | n        | Field contains (or does not contain) values
`starts_with`        | `not_starts_with`        | 1        | Field starts with (or does not start with) value
`ends_with`          | `not_ends_with`          | 1        | Field ends with (or does not end with) value
`between`            | `not_between`            | 2        | Field is between (or not between) values
`greater_than`       | `less_than`              | 1        | Field is greater than (or less than) value

Complex criteria can be constructed using `and`/`or` and parentheses.

To find derived records, the "->" operator can be used between two criteria. This will find all records that match the first criteria and then find all records that match the second criteria and are derived from the first (and so on). This operator can not be used inside parentheses. In order to limit the search space, parent records are expected to have the same ID as the derived records. This behaviour can be disabled using the `slims.unrestrict_parents` option. The `slims.bioinfo_check_criteria` currently does not support this operator.

eg. `cntn_foo equals a -> cntn_cstm_foo equals b and (cntn_bar not_between c d or cntn_baz equals e)`

## Fields

Fields are specified using the full SLIMS field name (eg. `cntn_fk_contentType`). Sub-fields inside JSON fields can be accessed by using the `json:` prefix and dot-notation (eg. json:cntn_cstm_foo.bar.baz)

## Map values

SLIMS fields can be mapped to keys on the `data.Samples` object using `<KEY>=<FIELD>`. This can be used to override any field on the `data.Samples` object (eg. `files=json:cntn_myJSONField.fastq_paths`).

> **NOTE** The key needs to be a valid attribute name on the `data.Samples` object. The `meta` attribute can be used to store arbitrary data (eg. `meta.my_key=cntn_SomeField`).

## Derivations

Derived records can be automatically created by specifying a mapping using `<FIELD>=<VALUE>` syntax (eg. `cntn_fk_contentType=1234 cntn_cstm_foo={state}`). Curly braces can be used to reference keys on the Sample object. The value of the field will be used as the value for the derived record. The derived record will be created if it does not exist and updated if it does. The derived record will be linked to the parent record using the `cntn_fk_originalContent` field. The `cntn_fk_contentType` field must be specified on the derived record. SLIMS may fail to create a records silently if required fields are not specified, and there is currently no way to explicitly check for required fields.

To create multiple derived records, multiple mappings can be specified as a list (or by passing the `--slims_derive` flag multiple times).