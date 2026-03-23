# SLIMS module for Cellophane

Fetch sample information from a SLIMS database. If samples are present (from a previous module or a samples file) these samples will be augmented with information from slims. If one or more sample IDs are specified they will be searched for. If no samples are present and no sample IDs are specified all samples matching the specified criteria will be fetched.

## Configuration

Option                            | Type      | Required | Default | Description
----------------------------------|-----------|----------|---------|-------------
`slims.url`                       | string    | Yes      |         | SLIMS URL
`slims.username`                  | string    | Yes      |         | SLIMS username
`slims.password`                  | string    | Yes      |         | SLIMS password
`slims.map`                       | mapping   |          |         | Mapping of keys to SLIMS field(s) (Use json: prefix and dot notation for JSON fields)
`slims.sync`                      | array     |          |         | Fields in slims.map that will be synced to SLIMS in a pre and post hook
`slims.match`                     | array     |          |         | Use fields from 'slims.map' when matching samples to augment with SLIMS metadata (in addition to 'sample.id')
`slims.derive`                    | mapping   |          |         | Mapping for creating derived records in SLIMS (Use curly braces to access keys in the 'sample' object, eg. '{sample.id}')
`slims.novel.max_age`             | string    |          | 1 year  | Maximum age of novel records to consider for matching (eg. "4 days", "1 year")
`slims.novel.criteria`            | string    |          |         | SLIMS criteria to use for matching novel records (In conjunction with 'slims.criteria')
`slims.criteria`                  | string    |          |         | SLIMS criteria to use for finding samples (eg. "cntn_cstm_SecondaryAnalysis equals 1337")
`slims.dry_run`                   | boolean   |          | false   | Do not sync data to SLIMS

### Example

```yaml
slims:
  # Ask your IT people for credentials (Remember to say please)
  url: 'http://slims.example.com/slimsrest'
  username: apiusername
  password: apipassword
  # Find sample records with content type 22 (Fastq) and qdrna in the samplesheet decription field
  # Exclude samples explicitly marked as "Do not include"
  criteria: |
    cntn_fk_contentType equals 22
    and cntn_cstm_doNotInclude not_equals true
    and cntn_cstm_rawSheetDescription contains qdrna
  novel:
    # Only consider samples that are less than 4 days old
    max_age: 4 days
    # Samples are novel if they do not have records with contentType 23 (Bioinformatics)
    criteria: |
      not_has_derived cntn_fk_contentType equals 23
  # Get file path from slims (unless specified in samples file)
  # Add backup remote keys for HCP hook
  # Add run information to samples
  map:
    files: json:cntn_cstm_demuxerSampleResult.fastq_paths
    run: cntn_cstm_runTag
    backup: json:cntn_cstm_demuxerBackupSampleResult.remote_keys
  # Include run tag when matching samples to SLIMS records
  match:
  - run
  # Create derived records with content type = 23 (Bioinformatics)
  # Include state from cellophane sample (access key using curly braces)
  # Other parameters are required by SLIMS or records will not be created
  derive:
    bioinformatics:
      cntn_fk_contentType: 23
      cntn_fk_location: 83
      cntn_cstm_SecondaryAnalysisState: '{sample.state}'
      cntn_status: 10
  sync:
  -

```

## Hooks

Name                   | When | Condition | Description
-----------------------|------|-----------|-------------
`slims_fetch`          | Pre  |           | Fetch sample info
`slims_sync_pre`       | Pre  |           | Sync data to SLIMS (pre hook)
`slims_sync_post`      | Post | Always    | Sync data to SLIMS (post hook)

> ![NOTE]
> The `slims_sync_pre` and `slims_sync_post` will sync any field specified in `slims.sync` to the original record based on `slims.map`. Any records specified in `slims.derive` will be created by the `slims_sync_pre` hook and updated by the `slims_sync_post` hook.

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
  criteria: str,
  config: cfg.Config,
  connection: slims.SLIMS | None = None,
) -> data.Samples
```

Classmethod to create a new data.Samples object with records matching the specified [criteria](criteria). Optionally pass a SLIMS connection object to avoid creating a new one.

```python
Samples.sync_derived(
  config: cfg.Config,
) -> None
```

Update bioinformatics state for all samples in SLIMS.

```python
Samples.sync_records(
  config: cfg.Config,
) -> None
```

Update the main record for all samples in SLIMS.

---

### Sample

```python
Sample.record: slims.Record | None
```

Record object for this sample (or None if sample has no Record).

```python
Sample.matches_record(
  record: Record,
  map_: dict,
  match: list[str] | None = None,
)
```

Check if the sample matches the specified record based on the `match` fields in the `map_` dictionary. Optionally pass a list of fields to match on using the `match` parameter.

```python
sample.match_from_record(
  record: slims.Record,
  map_: dict,
  map_ignore: list[str] | None = None,
```

Map values from a SLIMS record to the sample object based on the `map_` dictionary. Optionally pass a list of fields to ignore using the `map_ignore` parameter.

```python
Sample.from_record(
  record: slims.Record,
  config: cellophane.Config,
  **kwargs: Any,
)
```

Create a new data.Sample object from a single SLIMS record. Optionally pass additional keyword arguments to override fields on the sample object.

```python
Sample.sync_record(
  config: cellophane.Config,
)
```

Sync fields mapped to this sample object back to the SLIMS record.

```python
Sample.sync_derived(
  config: cellophane.Config,
)
```

Sync derived records for this sample back to SLIMS.

## Criteria

Criteria are specified using `<FIELD> <OPERATOR> <VALUE> ...` syntax where operators may take one or more values.

> [!WARNING]
> Invalid field names are silently ignored. This makes it very easy to accidentally fetch more samples than expected. This behavior is inherent to the official SLIMS API python bindings.

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
`has_parent`         | `not_has_parent`         | 1        | Field has (or does not have) a parent record that matches the specified criteria (where criteria is the value)
`has_derived`        | `not_has_derived`        | 1        | Field has (or does not have) a derived record that matches the specified criteria (where criteria is the value)

Complex criteria can be constructed using `and`/`or` and parentheses.

> ![IMPORTANT]
> The `has_parent` and `has_derived` operators are implemented as separate requests to the slims API. They should only be used to refine an already reasonably specific criteria. Unless the base criteria (without `has_parent`/`has_derived`) is specific enough this means a very large number of samples may be fetched, which can be very slow.

## Fields

Fields are specified using the full SLIMS field name (eg. `cntn_fk_contentType`). Sub-fields inside JSON fields can be accessed by using the `json:` prefix and dot-notation (eg. json:cntn_cstm_foo.bar.baz)

## Map values

SLIMS fields can be mapped to keys on the `data.Samples` object using `<KEY>=<FIELD>`. This can be used to override any field on the `data.Samples` object (eg. `files=json:cntn_myJSONField.fastq_paths`).

> ![NOTE]
> The key needs to be a valid attribute name on the `data.Samples` object. The `meta` attribute can be used to store arbitrary data (eg. `meta.my_key=cntn_SomeField`).

## Derivations

Derived records can be automatically created by specifying a mapping using `<FIELD>=<VALUE>` syntax (eg. `cntn_fk_contentType=1234 cntn_cstm_foo={state}`). Curly braces can be used to reference keys on the Sample object. The value of the field will be used as the value for the derived record. The derived record will be created if it does not exist and updated if it does. The derived record will be linked to the parent record using the `cntn_fk_originalContent` field. The `cntn_fk_contentType` field must be specified on the derived record. SLIMS may fail to create a records silently if required fields are not specified, and there is currently no way to explicitly check for required fields.

To create multiple derived records, multiple mappings can be specified as a list (or by passing the `--slims_derive` flag multiple times).