use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum PavanError {
    #[error("invalid angle: {0}")]
    InvalidAngle(String),
    #[error("invalid altitude: {0}")]
    InvalidAltitude(String),
    #[error("invalid velocity: {0}")]
    InvalidVelocity(String),
    #[error("invalid geometry: {0}")]
    InvalidGeometry(String),
    #[error("computation error: {0}")]
    ComputationError(String),
}

pub type Result<T> = std::result::Result<T, PavanError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_display() {
        let e = PavanError::InvalidAltitude("too high".into());
        assert!(e.to_string().contains("too high"));
    }

    #[test]
    fn result_type() {
        let ok: Result<f64> = Ok(1.0);
        assert!(ok.is_ok());
    }

    #[test]
    fn all_variants_display() {
        let cases: Vec<PavanError> = vec![
            PavanError::InvalidAngle("bad angle".into()),
            PavanError::InvalidAltitude("too high".into()),
            PavanError::InvalidVelocity("negative".into()),
            PavanError::InvalidGeometry("zero area".into()),
            PavanError::ComputationError("diverged".into()),
        ];
        let prefixes = [
            "invalid angle",
            "invalid altitude",
            "invalid velocity",
            "invalid geometry",
            "computation error",
        ];
        for (e, prefix) in cases.iter().zip(prefixes.iter()) {
            let msg = e.to_string();
            assert!(
                msg.starts_with(prefix),
                "expected prefix '{prefix}', got '{msg}'"
            );
        }
    }

    #[test]
    fn error_is_send_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PavanError>();
    }

    #[test]
    fn result_err_variant() {
        let err: Result<f64> = Err(PavanError::ComputationError("nan".into()));
        assert!(err.is_err());
        let Err(e) = err else { unreachable!() };
        assert!(matches!(e, PavanError::ComputationError(_)));
    }
}
