# Start/end mail module for Cellophane

Module for sending start/end mail notifications.

## Configuration

Option                 | Type | Required | Default             | Description
-----------------------|------|----------|---------------------|-------------
`mail.skip`            | bool |          | false               | Don't send mail
`mail.start.subject`   | str  |          | See [Start](#Start) | Subject of the mail (jinja2 template)
`mail.start.body`      | str  |          | See [Start](#Start) | Body of the mail (jinja2 template)
`mail.start.from_addr` | str  | -        |                     | From address to use
`mail.start.to_addr`   | list | -        |                     | List of recipients
`mail.start.cc_addr`   | list | -        |                     | List of CC recipients
`mail.end.subject`     | str  | -        | See [End](#End)     | Subject of the mail (jinja2 template)
`mail.end.body`        | str  | -        | See [End](#End)     | Body of the mail (jinja2 template)
`mail.smtp.host`       | str  | x        |                     | SMTP host
`mail.smtp.port`       | int  |          | 25                  | SMTP port
`mail.smtp.tls`        | bool |          | false               | Use TLS
`mail.smtp.user`       | str  |          |                     | SMTP username
`mail.smtp.password`   | str  |          |                     | SMTP password

## Hooks

Name                   | When | Condition | Description
-----------------------|------|-----------|-------------
`start_mail`           | Pre  |           | Send email when pipeline starts
`end_mail`             | Post | Always    | Send email when pipeline ends

## Start

**Subject:** `{{ analysis }} started`

**Body:**

```
{{ analysis }} has started for {{ samples.unique_ids|length }} sample(s).
The following samples are being analyzed:
{% for id in samples.unique_ids %}
{{ id }}
{%- endfor %}
```

## End

**Subject:** `{{ analysis }} finished`

**Body:**

```
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