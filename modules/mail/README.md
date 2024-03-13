# Start/end mail module for Cellophane

Module for sending start/end mail notifications.

## Configuration

Option                 | Type | Required | Default                 | Description
-----------------------|------|----------|-------------------------|-------------
mail.send              | bool |          | false                   | Send mail
mail.from_addr         | str  |          |                         | Default from address
mail.to_addr           | arr  |          |                         | Default list of recipients
mail.cc_addr           | arr  |          |                         | Default list of CC recipients
mail.start.subject     | str  |          | [Start Subject](#start) | Subject of the mail (jinja2 template)
mail.start.body        | str  |          | [Start Body](#start)    | Body of the mail (jinja2 template)
mail.end.subject       | str  |          | [End Subject](#end)     | Subject of the mail (jinja2 template)
mail.end.body          | str  |          | [End Body](#end)        | Body of the mail (jinja2 template)
mail.smtp.host         | str  | x        |                         | SMTP host
mail.smtp.port         | int  |          | 25                      | SMTP port
mail.smtp.tls          | bool |          | false                   | Use TLS
mail.smtp.user         | str  |          |                         | SMTP username
mail.smtp.password     | str  |          |                         | SMTP password

## Hooks

Name                   | When | Condition | Description
-----------------------|------|-----------|-------------
`start_mail`           | Pre  |           | Send email when pipeline starts
`end_mail`             | Post | Always    | Send email when pipeline ends

## Start

**Subject:** `{{ analysis }} started`

**Body:**

```jinja
{{ analysis }} has started for {{ samples.unique_ids|length }} sample(s).

The following samples are being analyzed:
{% for id in samples.unique_ids %}
{{ id }}
{%- endfor %}
```

## End

**Subject:** `{{ analysis }} finished`

**Body:**

```jinja
{{ analysis }} has finished processing {{ samples.unique_ids|length }} sample(s).
{% if samples.failed|length > 0 %}
❗️ Analysis failed for the following samples:
{% for id in samples.failed.unique_ids %}
{{ id }}
{%- endfor %}
{% endif %}
{%- if samples.complete|length > 0 %}
✅ Analysis completed successfully for the following samples:
{% for id in samples.complete.unique_ids %}
{{ id }}
{%- endfor %}
{% endif %}
```