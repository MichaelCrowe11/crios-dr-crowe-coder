# CriOS Discovery Engine - Deployment Guide

Complete deployment guide for the CriOS Discovery Engine with Nova UI integration.

## Quick Start (Development)

### Prerequisites
- Python 3.11+
- Node.js 18+
- Git

### 1. Clone and Setup Backend
```bash
git clone <repository>
cd CriOS

# Install Python dependencies
pip install -r requirements.txt

# Install additional dev dependencies
pip install -e ".[dev]"

# Run initial setup
make install-dev
```

### 2. Start Backend API
```bash
# Method 1: Using Make
make web

# Method 2: Direct uvicorn
python -m uvicorn crios.web.main:app --reload --port 8000

# Method 3: Using CLI
python -m crios.cli.main --help
```

### 3. Setup and Start Frontend
```bash
cd ui

# Install dependencies
npm install

# Start development server
npm run dev
```

### 4. Access Application
- **Frontend**: http://localhost:3000
- **API Documentation**: http://localhost:8000/docs
- **API Health**: http://localhost:8000/health

## Docker Deployment

### Single Command Setup
```bash
# Build and start all services
docker-compose up --build

# Start in background
docker-compose up -d --build
```

### Services
- **CriOS API**: http://localhost:8000
- **Nova UI**: http://localhost:3000
- **PostgreSQL**: localhost:5432
- **Redis**: localhost:6379

### Docker Architecture
```
┌─────────────────┐    ┌─────────────────┐
│   Nova UI       │    │   CriOS API     │
│   (Port 3000)   │────│   (Port 8000)   │
└─────────────────┘    └─────────────────┘
                              │
                    ┌─────────┴─────────┐
                    │                   │
            ┌───────────────┐   ┌───────────────┐
            │  PostgreSQL   │   │     Redis     │
            │  (Port 5432)  │   │  (Port 6379)  │
            └───────────────┘   └───────────────┘
```

## Production Deployment

### Environment Variables
```bash
# Backend (.env)
CRIOS_ENV=production
DATABASE_URL=postgresql://user:pass@host:5432/crios
REDIS_URL=redis://host:6379
SECRET_KEY=your-secret-key

# Frontend (.env.local)  
CRIOS_API_URL=https://api.your-domain.com
NODE_ENV=production
```

### SSL/HTTPS Setup
```bash
# Using Let's Encrypt with Certbot
certbot --nginx -d your-domain.com -d api.your-domain.com

# Update nginx config for SSL termination
# See nginx.conf example below
```

### Scaling Configuration
```yaml
# docker-compose.prod.yml
version: '3.8'
services:
  crios-api:
    deploy:
      replicas: 3
      resources:
        limits:
          cpus: '1.0'
          memory: 2G
        reservations:
          cpus: '0.5'
          memory: 1G
```

## API Integration

### Health Check Endpoints
```bash
# Backend health
curl http://localhost:8000/health

# Frontend health  
curl http://localhost:3000/api/health
```

### Available API Endpoints
```bash
# Molecular validation
POST /validate
Content-Type: application/json
{
  "smiles": "CCO",
  "standardize": true
}

# Crowe scoring
POST /crowe-score  
Content-Type: application/json
{
  "smiles": "CCO",
  "compound_origin": "synthetic",
  "therapeutic_areas": ["neurotherapeutic"]
}

# Ethics compliance
POST /ethics-check
Content-Type: application/json  
{
  "smiles": "CCO"
}

# Molecular similarity
POST /similarity
Content-Type: application/json
{
  "query_smiles": "CCO",
  "database_smiles": ["CCC", "CCCO"],
  "threshold": 0.7,
  "max_results": 10
}
```

## Development Workflow

### Backend Development
```bash
# Run tests
make test

# Run with hot reload
make web

# Code formatting
make format

# Type checking
make type-check

# Full quality checks
make qa
```

### Frontend Development
```bash
cd ui

# Development server with hot reload
npm run dev

# Type checking
npm run type-check

# Production build
npm run build

# Start production server
npm start
```

### Demo and Testing
```bash
# Run full system demo
make demo

# Run Python demo script
make demo-python

# Clean demo artifacts
make demo-clean
```

## Database Setup

### PostgreSQL (Optional)
```sql
-- Create database
CREATE DATABASE crios_discovery;
CREATE USER crios WITH ENCRYPTED PASSWORD 'your_password';
GRANT ALL PRIVILEGES ON DATABASE crios_discovery TO crios;

-- Initialize tables
\i sql/schema.sql
```

### Redis (Optional)
```bash
# Start Redis
redis-server

# Or with Docker
docker run -d -p 6379:6379 redis:7-alpine
```

## Monitoring and Logging

### Application Logs
```bash
# Backend logs
tail -f logs/crios-api.log

# Frontend logs
docker logs crios-frontend

# API access logs
tail -f logs/access.log
```

### Health Monitoring
```bash
# API health check
curl -f http://localhost:8000/health || echo "API DOWN"

# Frontend health check  
curl -f http://localhost:3000 || echo "UI DOWN"
```

## Troubleshooting

### Common Issues

1. **RDKit Installation Issues**
   ```bash
   # On macOS with conda
   conda install -c conda-forge rdkit
   
   # On Ubuntu/Debian
   pip install rdkit --extra-index-url https://pypi.anaconda.org/conda-forge/simple
   ```

2. **Port Conflicts**
   ```bash
   # Check what's using port 8000
   lsof -i :8000
   
   # Kill process if needed
   kill -9 <PID>
   ```

3. **API Connection Issues**
   ```bash
   # Check API is running
   curl http://localhost:8000/health
   
   # Check frontend environment
   echo $CRIOS_API_URL
   ```

4. **Docker Issues**
   ```bash
   # Rebuild containers
   docker-compose down && docker-compose up --build
   
   # Clear Docker cache
   docker system prune -a
   ```

### Performance Optimization

1. **Backend Optimization**
   ```bash
   # Enable multiple workers
   gunicorn crios.web.main:app -w 4 -k uvicorn.workers.UvicornWorker
   ```

2. **Frontend Optimization**
   ```bash
   # Enable compression
   npm run build
   
   # Analyze bundle size
   npm install --save-dev @next/bundle-analyzer
   ```

## Security Considerations

### API Security
- Enable CORS properly for production
- Implement rate limiting
- Use HTTPS in production
- Validate all inputs
- Implement authentication if needed

### Frontend Security
- Sanitize user inputs
- Use environment variables for sensitive config
- Implement CSP headers
- Regular dependency updates

## Backup and Recovery

### Data Backup
```bash
# PostgreSQL backup
pg_dump crios_discovery > backup_$(date +%Y%m%d).sql

# Redis backup
redis-cli BGSAVE
```

### Configuration Backup
```bash
# Backup configuration files
tar -czf config_backup_$(date +%Y%m%d).tar.gz configs/ .env ui/.env.local
```

## Fly.io Cloud Deployment

### Prerequisites
```bash
# Install Fly.io CLI
curl -L https://fly.io/install.sh | sh

# Login to Fly.io
flyctl auth login
```

### 1. Create PostgreSQL Database
```bash
flyctl postgres create --name crios-db --region iad --vm-size shared-cpu-1x --volume-size 10
```

### 2. Create Redis Instance
```bash
flyctl redis create --name crios-redis --region iad --plan free
```

### 3. Deploy Backend
```bash
cd platform/backend

# Create the app (first time only)
flyctl apps create crios-backend --org personal

# Set secrets
flyctl secrets set \
  DATABASE_URL="postgres://username:password@hostname:5432/database" \
  REDIS_URL="redis://default:password@hostname:6379" \
  ANTHROPIC_API_KEY="your-api-key"

# Deploy
flyctl deploy

# Check status
flyctl status
```

### 4. Deploy Frontend
```bash
cd platform/frontend

# Create the app (first time only)
flyctl apps create crios-frontend --org personal

# Set environment variables
flyctl secrets set \
  NEXT_PUBLIC_API_URL="https://crios-backend.fly.dev" \
  NEXT_PUBLIC_WS_URL="wss://crios-backend.fly.dev"

# Deploy
flyctl deploy
```

### Fly.io Configuration Files
- Backend: `platform/backend/fly.toml`
- Frontend: `platform/frontend/fly.toml`

### Monitoring on Fly.io
```bash
# View logs
flyctl logs -a crios-backend
flyctl logs -a crios-frontend

# Check status
flyctl status -a crios-backend

# SSH into machine
flyctl ssh console -a crios-backend

# Scale resources
flyctl scale count 2 -a crios-backend
flyctl scale memory 2048 -a crios-backend
```

### URLs After Fly.io Deployment
- **Backend API**: https://crios-backend.fly.dev
- **Frontend**: https://crios-frontend.fly.dev
- **API Health**: https://crios-backend.fly.dev/health

### Troubleshooting Fly.io
1. **Deployment fails**: Check logs with `flyctl logs`
2. **Out of memory**: Scale up with `flyctl scale memory 2048`
3. **Connection issues**: Verify secrets with `flyctl secrets list`
4. **Rollback**: Use `flyctl releases rollback <version>`

## Support

- **Documentation**: See README.md files in each component
- **API Docs**: http://localhost:8000/docs
- **Issues**: Create GitHub issues for bugs/features
- **Development**: Use `make help` for available commands
- **Fly.io Docs**: https://fly.io/docs/

---

**CriOS Discovery Engine** • Science before status. Discovery before profit. • © 2025 Crowe BioSystems