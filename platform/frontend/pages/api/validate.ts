import type { NextApiRequest, NextApiResponse } from 'next';

type ValidateRequest = {
  smiles: string;
  standardize?: boolean;
};

type ValidateResponse = {
  valid: boolean;
  smiles?: string;
  standardized_smiles?: string;
  molecular_weight?: number;
  error?: string;
};

/**
 * Molecule validation endpoint
 * POST /api/validate
 *
 * Note: This is a lightweight version for Vercel.
 * For full RDKit functionality, consider using the Python API route or external service.
 */
export default async function handler(
  req: NextApiRequest,
  res: NextApiResponse<ValidateResponse>
) {
  // Only allow POST requests
  if (req.method !== 'POST') {
    res.setHeader('Allow', ['POST']);
    return res.status(405).json({
      valid: false,
      error: 'Method not allowed. Use POST.',
    });
  }

  try {
    const { smiles, standardize = false }: ValidateRequest = req.body;

    if (!smiles) {
      return res.status(400).json({
        valid: false,
        error: 'SMILES string is required',
      });
    }

    // Basic SMILES validation (simplified)
    // In production, this should call RDKit or external chemistry API
    const isValid = validateSmilesBasic(smiles);

    if (!isValid) {
      return res.status(400).json({
        valid: false,
        smiles,
        error: 'Invalid SMILES string',
      });
    }

    // For actual molecular validation, you would:
    // 1. Use a Python serverless function with RDKit
    // 2. Call an external chemistry API
    // 3. Use a lightweight JavaScript chemistry library

    res.status(200).json({
      valid: true,
      smiles,
      standardized_smiles: standardize ? smiles : undefined,
      molecular_weight: undefined, // Would be calculated with RDKit
    });

  } catch (error) {
    console.error('Validation error:', error);
    res.status(500).json({
      valid: false,
      error: error instanceof Error ? error.message : 'Internal server error',
    });
  }
}

/**
 * Basic SMILES validation
 * This is a simplified version - use RDKit for production
 */
function validateSmilesBasic(smiles: string): boolean {
  if (!smiles || typeof smiles !== 'string') {
    return false;
  }

  // Basic checks
  const validChars = /^[A-Za-z0-9@+\-\[\]()=#$:\/\\.*%]+$/;
  if (!validChars.test(smiles)) {
    return false;
  }

  // Check for balanced parentheses
  let parenCount = 0;
  let bracketCount = 0;

  for (const char of smiles) {
    if (char === '(') parenCount++;
    if (char === ')') parenCount--;
    if (char === '[') bracketCount++;
    if (char === ']') bracketCount--;

    if (parenCount < 0 || bracketCount < 0) {
      return false;
    }
  }

  return parenCount === 0 && bracketCount === 0;
}
