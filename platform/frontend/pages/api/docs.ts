import type { NextApiRequest, NextApiResponse } from 'next';

type ApiDocsResponse = {
  title: string;
  version: string;
  description: string;
  endpoints: {
    [key: string]: {
      method: string;
      path: string;
      description: string;
      parameters?: any;
      responses?: any;
    };
  };
};

/**
 * API documentation endpoint
 * GET /api/docs
 */
export default async function handler(
  req: NextApiRequest,
  res: NextApiResponse<ApiDocsResponse | string>
) {
  if (req.method !== 'GET') {
    res.setHeader('Allow', ['GET']);
    return res.status(405).end('Method Not Allowed');
  }

  // Return API documentation
  const docs: ApiDocsResponse = {
    title: 'CriOS Discovery Engine API',
    version: '1.0.0',
    description: 'AI-Driven Drug Discovery Platform - Vercel Deployment',
    endpoints: {
      health: {
        method: 'GET',
        path: '/api/health',
        description: 'Health check endpoint',
        responses: {
          200: {
            description: 'Service is healthy',
            schema: {
              status: 'string',
              timestamp: 'string',
              service: 'string',
              version: 'string',
            },
          },
        },
      },
      validate: {
        method: 'POST',
        path: '/api/validate',
        description: 'Validate molecular structure (SMILES)',
        parameters: {
          smiles: 'string (required) - SMILES representation of molecule',
          standardize: 'boolean (optional) - Whether to standardize the molecule',
        },
        responses: {
          200: {
            description: 'Validation successful',
            schema: {
              valid: 'boolean',
              smiles: 'string',
              standardized_smiles: 'string',
              molecular_weight: 'number',
            },
          },
          400: {
            description: 'Invalid input',
          },
        },
      },
      discovery: {
        method: 'POST',
        path: '/api/discovery/start',
        description: 'Start drug discovery pipeline',
        parameters: {
          target: 'string - Disease target',
          compounds: 'array - Initial compound library',
          pipeline: 'string - Pipeline type (full, fast, custom)',
        },
      },
      agents: {
        method: 'GET',
        path: '/api/agents/status',
        description: 'Get status of all PhD agents',
      },
    },
  };

  // Check if HTML documentation is requested
  const acceptHeader = req.headers.accept || '';
  if (acceptHeader.includes('text/html')) {
    // Return HTML version
    const html = generateHtmlDocs(docs);
    res.setHeader('Content-Type', 'text/html');
    return res.status(200).send(html);
  }

  // Return JSON by default
  res.status(200).json(docs);
}

function generateHtmlDocs(docs: ApiDocsResponse): string {
  return `
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>${docs.title} - API Documentation</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      max-width: 1200px;
      margin: 0 auto;
      padding: 20px;
      background: #0a0a0a;
      color: #e0e0e0;
    }
    h1 { color: #00ff88; border-bottom: 2px solid #00ff88; padding-bottom: 10px; }
    h2 { color: #00ccff; margin-top: 30px; }
    .endpoint {
      background: #1a1a1a;
      border-left: 4px solid #00ff88;
      padding: 15px;
      margin: 20px 0;
      border-radius: 4px;
    }
    .method {
      display: inline-block;
      padding: 4px 12px;
      border-radius: 4px;
      font-weight: bold;
      font-size: 12px;
      margin-right: 10px;
    }
    .get { background: #00ff88; color: #000; }
    .post { background: #00ccff; color: #000; }
    .path { font-family: 'Courier New', monospace; color: #ffcc00; }
    code { background: #2a2a2a; padding: 2px 6px; border-radius: 3px; color: #00ff88; }
    pre { background: #2a2a2a; padding: 15px; border-radius: 4px; overflow-x: auto; }
  </style>
</head>
<body>
  <h1>${docs.title}</h1>
  <p><strong>Version:</strong> ${docs.version}</p>
  <p>${docs.description}</p>

  <h2>Endpoints</h2>
  ${Object.entries(docs.endpoints)
    .map(
      ([key, endpoint]) => `
    <div class="endpoint">
      <div>
        <span class="method ${endpoint.method.toLowerCase()}">${endpoint.method}</span>
        <span class="path">${endpoint.path}</span>
      </div>
      <p>${endpoint.description}</p>
      ${
        endpoint.parameters
          ? `
        <h3>Parameters</h3>
        <pre>${JSON.stringify(endpoint.parameters, null, 2)}</pre>
      `
          : ''
      }
      ${
        endpoint.responses
          ? `
        <h3>Responses</h3>
        <pre>${JSON.stringify(endpoint.responses, null, 2)}</pre>
      `
          : ''
      }
    </div>
  `
    )
    .join('')}

  <footer style="margin-top: 50px; padding-top: 20px; border-top: 1px solid #333; color: #666;">
    <p>© 2025 CriOS Discovery Engine | Science before status. Discovery before profit.</p>
  </footer>
</body>
</html>
  `.trim();
}
