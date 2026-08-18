use std::collections::HashMap;

use csv::StringRecord;

use crate::target_list::{TargetErrorInner, TargetErrors};

pub(super) fn read_csv_trimmed(target_list: &str) -> csv::Reader<&[u8]> {
    let target_list = target_list.trim();

    csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(target_list.as_bytes())
}

pub(super) fn rename_fields(
    original_fieldnames: &StringRecord,
    field_aliases: &HashMap<&str, &str>,
) -> (StringRecord, Option<TargetErrors>) {
    let mut renamed_fields = StringRecord::new();
    let mut errors = Vec::new();

    for original in original_fieldnames {
        let renamed = field_aliases.get(original).unwrap_or(&original);

        renamed_fields.push_field(renamed);

        if *renamed != original {
            errors.push(TargetErrorInner::RenamedField {
                original_fieldname: original.to_owned(),
                correct_fieldname: (*renamed).to_owned(),
            });
        }
    }

    (
        renamed_fields,
        (!errors.is_empty()).then_some(TargetErrors {
            line_number: None,
            submitted_target: None,
            errors: errors.into_iter().map(Into::into).collect(),
        }),
    )
}

#[cfg(test)]
mod tests {
    use crate::target_list::{TargetErrorInner, TargetErrors, csv_util::rename_fields};

    #[test]
    fn renaming_fields() {
        let original_fieldnames = ["field1", "field2"].iter().collect();
        let field_aliases = [("field1", "field_1")].into_iter().collect();

        let (renamed_fields, error) = rename_fields(&original_fieldnames, &field_aliases);

        assert_eq!(
            renamed_fields,
            ["field_1", "field2"][..],
            "failed to rename fields"
        );

        assert_eq!(
            error,
            Some(TargetErrors {
                line_number: None,
                submitted_target: None,
                errors: vec![
                    TargetErrorInner::RenamedField {
                        original_fieldname: "field1".to_owned(),
                        correct_fieldname: "field_1".to_owned()
                    }
                    .into()
                ]
            }),
            "failed to construct field-renaming error"
        );
    }
}
