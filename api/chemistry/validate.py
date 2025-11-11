"""
Molecule validation endpoint using Python serverless function
Note: RDKit is too large for Vercel serverless (>250MB limit)

For production, consider:
1. Using external chemistry API (PubChem, ChemSpider)
2. Deploying RDKit backend separately (Railway, Render)
3. Using CDK (Chemistry Development Kit) via Java/Node
4. Custom lightweight SMILES parser
"""

from http.server import BaseHTTPRequestHandler
import json
from typing import Dict, Any


class handler(BaseHTTPRequestHandler):
    """
    Serverless function handler for Vercel
    POST /api/chemistry/validate
    """

    def do_POST(self):
        try:
            # Read request body
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length)
            data = json.loads(body.decode('utf-8'))

            smiles = data.get('smiles', '')
            standardize = data.get('standardize', False)

            if not smiles:
                self.send_error_response(400, 'SMILES string is required')
                return

            # Validate SMILES (basic validation without RDKit)
            result = validate_smiles_lightweight(smiles, standardize)

            self.send_success_response(result)

        except json.JSONDecodeError:
            self.send_error_response(400, 'Invalid JSON')
        except Exception as e:
            self.send_error_response(500, str(e))

    def send_success_response(self, data: Dict[str, Any]):
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        self.wfile.write(json.dumps(data).encode('utf-8'))

    def send_error_response(self, code: int, message: str):
        self.send_response(code)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        error_data = {'error': message, 'valid': False}
        self.wfile.write(json.dumps(error_data).encode('utf-8'))


def validate_smiles_lightweight(smiles: str, standardize: bool = False) -> Dict[str, Any]:
    """
    Lightweight SMILES validation without RDKit
    For production, replace with actual RDKit validation or external API call
    """

    # Basic validation
    if not isinstance(smiles, str) or not smiles:
        return {'valid': False, 'error': 'Invalid SMILES string'}

    # Check for valid characters
    valid_chars = set('ABCNOPSFIKabcnopsfi0123456789@+\\-[]()=#$:/%.*')
    if not all(c in valid_chars for c in smiles):
        return {'valid': False, 'error': 'Invalid characters in SMILES'}

    # Check balanced brackets
    if not check_balanced_brackets(smiles):
        return {'valid': False, 'error': 'Unbalanced brackets in SMILES'}

    # Return validation result
    return {
        'valid': True,
        'smiles': smiles,
        'standardized_smiles': smiles if standardize else None,
        'note': 'Using lightweight validation. For production, use RDKit or external API.',
        'recommendation': 'Deploy RDKit backend on Railway or use PubChem API for full validation'
    }


def check_balanced_brackets(smiles: str) -> bool:
    """Check if brackets are balanced"""
    paren_count = 0
    bracket_count = 0

    for char in smiles:
        if char == '(':
            paren_count += 1
        elif char == ')':
            paren_count -= 1
        elif char == '[':
            bracket_count += 1
        elif char == ']':
            bracket_count -= 1

        if paren_count < 0 or bracket_count < 0:
            return False

    return paren_count == 0 and bracket_count == 0
