use thiserror::Error;

#[derive(Error, Debug)]
pub enum MolioError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Unsupported format: {0}")]
    UnsupportedFormat(String),

    #[error("Invalid atom data: {0}")]
    InvalidAtom(String),

    #[error("Invalid bond data: {0}")]
    InvalidBond(String),
}

pub type MolioResult<T> = Result<T, MolioError>;
