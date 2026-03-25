/// Initialize tracing subscriber for pavan logging.
pub fn init() {
    use tracing_subscriber::EnvFilter;
    let filter = EnvFilter::try_from_env("PAVAN_LOG").unwrap_or_else(|_| EnvFilter::new("warn"));
    tracing_subscriber::fmt().with_env_filter(filter).init();
}
