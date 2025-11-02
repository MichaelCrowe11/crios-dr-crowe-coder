# Railway Deployment Guide for CriOS Dr. Crowe Coder

Complete guide to deploying the CriOS Dr. Crowe Coder platform on Railway with all services integrated.

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Quick Start](#quick-start)
3. [Service Architecture](#service-architecture)
4. [Step-by-Step Deployment](#step-by-step-deployment)
5. [Environment Variables](#environment-variables)
6. [Database Setup](#database-setup)
7. [Monitoring & Debugging](#monitoring--debugging)
8. [Scaling & Optimization](#scaling--optimization)
9. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required
- **Railway Account**: Sign up at [railway.app](https://railway.app)
- **Railway CLI**: Install via `npm i -g @railway/cli` or `brew install railway`
- **Git**: Ensure your code is in a Git repository
- **Anthropic API Key**: Get from [console.anthropic.com](https://console.anthropic.com)

### Optional (for advanced features)
- Stripe account (for payments)
- SendGrid account (for emails)
- AWS S3 bucket (for file storage)
- Sentry account (for error monitoring)

---

## Quick Start

### 1. Install Railway CLI
```bash
# npm
npm i -g @railway/cli

# Homebrew
brew install railway

# Verify installation
railway --version
```

### 2. Login to Railway
```bash
railway login
```

### 3. Create a New Project
```bash
# From the project root
cd /path/to/crios-dr-crowe-coder

# Initialize Railway project
railway init

# Link to Railway
railway link
```

### 4. Deploy All Services
```bash
# Deploy using Railway CLI
railway up

# Or link your GitHub repository for automatic deployments
```

---

## Service Architecture

The CriOS platform consists of **5 main services**:

```
┌─────────────────────────────────────────────────────┐
│                  Railway Project                    │
├─────────────────────────────────────────────────────┤
│  1. PostgreSQL (Railway-managed)                    │
│  2. Redis (Railway-managed)                         │
│  3. Backend API (FastAPI + RDKit)                   │
│  4. Frontend (Next.js 14)                           │
│  5. Celery Worker (Background Tasks)                │
│  6. Core API (CriOS Discovery Engine)               │
└─────────────────────────────────────────────────────┘
```

### Service Details

| Service | Description | Port | Health Check |
|---------|-------------|------|--------------|
| **PostgreSQL** | Main database | 5432 | Managed by Railway |
| **Redis** | Cache & message broker | 6379 | Managed by Railway |
| **Backend API** | FastAPI server | 8000 | `/health` |
| **Frontend** | Next.js dashboard | 3000 | `/` |
| **Celery Worker** | ML pipeline & background jobs | N/A | N/A |
| **Core API** | CriOS Discovery Engine | 8000 | `/health` |

---

## Step-by-Step Deployment

### Step 1: Create Railway Project
```bash
# Login
railway login

# Create new project
railway init --name crios-dr-crowe-coder

# Or use the Railway dashboard:
# 1. Go to railway.app
# 2. Click "New Project"
# 3. Choose "Empty Project"
```

### Step 2: Add PostgreSQL Database
```bash
# Via CLI
railway add --plugin postgresql

# Or via Dashboard:
# 1. Click "+ New"
# 2. Select "Database"
# 3. Choose "PostgreSQL"
```

Railway will automatically:
- Provision a PostgreSQL instance
- Create environment variables (`DATABASE_URL`, `POSTGRES_USER`, etc.)
- Make them available to all services

### Step 3: Add Redis
```bash
# Via CLI
railway add --plugin redis

# Or via Dashboard:
# 1. Click "+ New"
# 2. Select "Database"
# 3. Choose "Redis"
```

### Step 4: Deploy Backend API
```bash
# Create backend service
railway service create crios-backend

# Set build configuration
railway service configure crios-backend \
  --dockerfile platform/backend/Dockerfile \
  --root platform/backend

# Add environment variables (see .env.railway.template)
railway variables set ANTHROPIC_API_KEY=your_key_here

# Deploy
railway up --service crios-backend
```

**Environment Variables for Backend:**
```bash
ANTHROPIC_API_KEY=sk-ant-...
DATABASE_URL=${{Postgres.DATABASE_URL}}
REDIS_URL=${{Redis.REDIS_URL}}
JWT_SECRET=your_long_random_secret
PYTHONPATH=/app
PYTHONUNBUFFERED=1
WORKERS=4
```

### Step 5: Deploy Frontend
```bash
# Create frontend service
railway service create crios-frontend

# Set build configuration
railway service configure crios-frontend \
  --dockerfile platform/frontend/Dockerfile \
  --root platform/frontend

# Add environment variables
railway variables set \
  NODE_ENV=production \
  NEXT_PUBLIC_API_URL=https://crios-backend.railway.app \
  NEXT_TELEMETRY_DISABLED=1

# Deploy
railway up --service crios-frontend
```

**Environment Variables for Frontend:**
```bash
NODE_ENV=production
NEXT_TELEMETRY_DISABLED=1
NEXT_PUBLIC_API_URL=https://${{crios-backend.RAILWAY_PUBLIC_DOMAIN}}
NEXT_PUBLIC_WS_URL=wss://${{crios-backend.RAILWAY_PUBLIC_DOMAIN}}
```

### Step 6: Deploy Celery Worker
```bash
# Create worker service
railway service create crios-worker

# Set build configuration
railway service configure crios-worker \
  --dockerfile platform/backend/Dockerfile.worker \
  --root platform/backend

# Use same env vars as backend, plus:
railway variables set C_FORCE_ROOT=true

# Deploy
railway up --service crios-worker
```

### Step 7: Deploy Core API (Optional)
```bash
# Create core API service
railway service create crios-core-api

# Set build configuration
railway service configure crios-core-api \
  --dockerfile Dockerfile.api \
  --root .

# Deploy
railway up --service crios-core-api
```

---

## Environment Variables

### Using Railway Dashboard
1. Go to your project
2. Select a service
3. Click "Variables" tab
4. Add variables from `.env.railway.template`

### Using Railway CLI
```bash
# Set single variable
railway variables set KEY=value

# Set multiple variables
railway variables set \
  KEY1=value1 \
  KEY2=value2 \
  KEY3=value3

# Set from file
railway variables set --file .env.railway

# View all variables
railway variables
```

### Critical Variables

**Required for all services:**
- `ANTHROPIC_API_KEY`: Your Claude API key
- `DATABASE_URL`: Auto-set by Railway when PostgreSQL is added
- `REDIS_URL`: Auto-set by Railway when Redis is added

**Backend-specific:**
- `JWT_SECRET`: Random 32+ character string for auth
- `PYTHONPATH=/app`
- `PYTHONUNBUFFERED=1`

**Frontend-specific:**
- `NEXT_PUBLIC_API_URL`: URL of your backend service
- `NODE_ENV=production`

---

## Database Setup

### Automatic Setup
Railway automatically creates and connects PostgreSQL when you add it.

### Manual Migration (if needed)
```bash
# SSH into backend service
railway run --service crios-backend bash

# Run migrations
python -m alembic upgrade head

# Or use the backend app's migration command
python app.py migrate
```

### Database Access
```bash
# Connect to PostgreSQL
railway connect postgres

# Or get connection string
railway variables get DATABASE_URL
```

---

## Monitoring & Debugging

### View Logs
```bash
# All services
railway logs

# Specific service
railway logs --service crios-backend

# Follow logs (live)
railway logs --service crios-backend --follow
```

### Check Service Status
```bash
# List all services
railway status

# Get service details
railway service info crios-backend
```

### SSH into Container
```bash
# Access running container
railway run --service crios-backend bash

# Run one-off command
railway run --service crios-backend python manage.py shell
```

### Health Checks
```bash
# Backend API
curl https://your-backend-url.railway.app/health

# Core API
curl https://your-core-api-url.railway.app/health

# Frontend
curl https://your-frontend-url.railway.app/
```

---

## Scaling & Optimization

### Vertical Scaling
Railway automatically scales based on your plan:
- **Hobby**: 512MB RAM, 1 vCPU
- **Pro**: 8GB RAM, 8 vCPU
- **Enterprise**: Custom resources

### Horizontal Scaling
```bash
# Increase replicas (Pro+ plan)
railway service scale crios-backend --replicas 3
```

### Performance Tips

1. **Enable Redis Caching**
   - Railway's Redis is optimized for production
   - Cache expensive ML predictions

2. **Optimize Docker Images**
   - Multi-stage builds (already configured)
   - Layer caching for faster deploys

3. **Database Indexing**
   - Add indexes for frequently queried fields
   - Use connection pooling

4. **CDN for Frontend**
   - Railway provides global CDN
   - Enable in service settings

---

## Troubleshooting

### Common Issues

#### 1. Build Failures
```bash
# Check build logs
railway logs --service crios-backend --deployment latest

# Common fixes:
# - Ensure Dockerfile paths are correct
# - Check requirements.txt has all dependencies
# - Verify system packages are installed
```

#### 2. RDKit Installation Issues
```bash
# RDKit requires system libraries
# Ensure Dockerfile includes:
RUN apt-get install -y \
    librdkit-dev \
    libboost-all-dev
```

#### 3. Database Connection Issues
```bash
# Verify DATABASE_URL is set
railway variables get DATABASE_URL

# Check PostgreSQL is running
railway status postgres

# Test connection
railway run --service crios-backend python -c "import psycopg2; print('OK')"
```

#### 4. Frontend Can't Connect to Backend
```bash
# Ensure NEXT_PUBLIC_API_URL is set correctly
railway variables get NEXT_PUBLIC_API_URL

# Should be: https://crios-backend.railway.app (or your custom domain)

# Update and redeploy
railway variables set NEXT_PUBLIC_API_URL=https://your-backend-url.railway.app
railway up --service crios-frontend
```

#### 5. Celery Worker Not Processing Jobs
```bash
# Check worker logs
railway logs --service crios-worker --follow

# Ensure REDIS_URL is set
railway variables get REDIS_URL

# Verify Redis is accessible
railway run --service crios-worker redis-cli -u $REDIS_URL ping
```

---

## Custom Domains

### Add Custom Domain
```bash
# Via CLI
railway domain add yourdomain.com --service crios-frontend

# Via Dashboard:
# 1. Go to service settings
# 2. Click "Domains"
# 3. Add custom domain
# 4. Update DNS records as instructed
```

### SSL Certificates
Railway automatically provisions SSL certificates for:
- Railway-provided domains (*.railway.app)
- Custom domains (via Let's Encrypt)

---

## Deployment Workflow

### Automatic Deployments (Recommended)
1. Connect GitHub repository in Railway dashboard
2. Choose branch (e.g., `main` or `production`)
3. Railway auto-deploys on every push
4. Set up branch-based environments

### Manual Deployments
```bash
# Deploy from local
railway up

# Deploy specific service
railway up --service crios-backend

# Deploy from specific branch
git checkout production
railway up
```

---

## Cost Optimization

### Railway Pricing
- **Hobby**: $5/month + $0.000231/GB-hour
- **Pro**: $20/month + $0.000231/GB-hour
- **Enterprise**: Custom pricing

### Tips to Reduce Costs
1. **Stop unused services**: Railway charges for running time
2. **Use shared databases**: One PostgreSQL/Redis for multiple services
3. **Optimize Docker images**: Smaller images = faster builds = lower costs
4. **Scale down development environments**: Use smaller instances for dev/staging

---

## Backup & Disaster Recovery

### Database Backups
```bash
# Railway provides automatic daily backups (Pro+ plan)
# Manual backup:
railway run --service crios-backend pg_dump $DATABASE_URL > backup.sql

# Restore from backup:
railway run --service crios-backend psql $DATABASE_URL < backup.sql
```

### Environment Variable Backups
```bash
# Export all variables
railway variables > env-backup.txt

# Restore variables
cat env-backup.txt | xargs railway variables set
```

---

## Next Steps

After successful deployment:

1. ✅ Test all endpoints
2. ✅ Run health checks
3. ✅ Set up monitoring (Sentry, Grafana)
4. ✅ Configure custom domains
5. ✅ Set up CI/CD pipeline
6. ✅ Enable automatic backups
7. ✅ Document API endpoints
8. ✅ Set up staging environment

---

## Support & Resources

- **Railway Docs**: https://docs.railway.app
- **Railway Discord**: https://discord.gg/railway
- **CriOS GitHub**: https://github.com/michaelcrowe11/crios
- **Railway Status**: https://status.railway.app

---

## Deployment Checklist

### Pre-Deployment
- [ ] Railway account created
- [ ] Railway CLI installed and logged in
- [ ] Anthropic API key obtained
- [ ] Environment variables prepared
- [ ] Code pushed to Git repository
- [ ] Dockerfiles tested locally

### Deployment
- [ ] PostgreSQL database provisioned
- [ ] Redis cache provisioned
- [ ] Backend API deployed and healthy
- [ ] Frontend deployed and accessible
- [ ] Celery worker running
- [ ] Environment variables set for all services
- [ ] Database migrations run

### Post-Deployment
- [ ] Health checks passing
- [ ] Frontend can communicate with backend
- [ ] Background jobs processing correctly
- [ ] Custom domains configured (if applicable)
- [ ] SSL certificates active
- [ ] Monitoring enabled
- [ ] Backups configured
- [ ] Documentation updated

---

**Deployment Time Estimate**: 30-45 minutes for complete setup

**Good luck with your deployment! 🚀**
