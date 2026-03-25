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
}
